#
# This file is part of the magnum.np distribution
# (https://gitlab.com/magnum.np/magnum.np).
# Copyright (c) 2023 magnum.np team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

from magnumnp.common import logging, timedmethod, constants, normalize
from . import LLGSolver
from magnumnp.linear_elasticity.bcs import Plane, PlaneBC
from magnumnp.linear_elasticity.deriv_term_compiler import *
from magnumnp.linear_elasticity.strain import epsilon, epsilon_m, _get_C_jump_conditions, _get_B_jump_conditions, _get_sigM_jump_conditions
from magnumnp.linear_elasticity.stress import sigma
from magnumnp.linear_elasticity.utils import gradient_with_pbc, _get_diff_slices

from .ode_solvers import RKF45
import torch
import numpy as np

__all__ = ["LLGWithLESolver"]

class LLGWithLESolver(LLGSolver):
    r"""
    Extension of the ``LLGSolver`` class that also yields the rate of change of 
    mechanical displacement and momentum density.

    When used with the magnetoelastic field, this solver enables self-consistent 
    time integration of magnetization dynamics and elastodynamics.

    Examples
    --------
    Standard usage:

        .. code:: python

            llg = LLGWithLESolver([magnetoelastic, demag, exchange, external])
            logger = Logger("data", fields=['m', 'ud', 'pd'])
            while state.t < 1e-9 - eps:
                llg.step(state, 1e-11)
                logger << state

    Optimizations for speed and lower memory usage:

    1. Cubic material with a magnetic film in the 3rd x-layer:

        .. code:: python

            llg = LLGWithLESolver(
                [magnetoelastic, anisotropy, demag, exchange, external],
                C_sym="cubic",
                magnetic_x_limits=[3, 4]
            )

    2. Material with a magnetic domain on [:, 2:, :4]:

        .. code:: python

            llg = LLGWithLESolver(
                [magnetoelastic, demag, exchange, external],
                magnetic_y_limits=[3, None],
                magnetic_z_limits=[0, 4],
                boundary_nodes=2
            )

    3. 1D setup with pbcs:

        .. code:: python

            llg = LLGWithLESolver(
                [magnetoelastic, anisotropy, exchange, external],
                iteration_depth = 0
            )

    Parameters
    ----------
    terms : list of :class:`LLGTerm`
        List of LLG contributions to be considered for time integration.
    solver : :class:`Solver`, optional
        ODE solver to use. Choose one of: ``RKF45`` (default), ``RKF56``, 
        ``ScipyODE``, ``ScipyOdeint``, ``TorchDiffEq``, or ``TorchDiffEqAdjoint``.
    no_precession : bool, optional
        If ``True``, integrates without the precession term (default: ``False``).
    mask_elastic : torch.Tensor, optional
        Boolean or integer mask of shape ``mesh.n`` restricting the region where 
        elastodynamic properties are updated (default: ``None``).
    C_sym : numpy.ndarray or str, optional
        Boolean or integer mask of shape ``(6, 6)`` restricting which stiffness 
        tensor elements are considered (default: ``None``). Can also be a string 
        (``"isotropic"`` or ``"cubic"``) to include only C11, C12, C13, C22, C23, 
        C33, C44, C55, and C66.
    magnetic_x_limits : list of int, optional
        Two-element list giving start and stop indices for slicing the x-direction 
        of the magnetic domain (default: ``None``).
    magnetic_y_limits : list of int, optional
        Two-element list giving start and stop indices for slicing the y-direction 
        of the magnetic domain (default: ``None``).
    magnetic_z_limits : list of int, optional
        Two-element list giving start and stop indices for slicing the z-direction 
        of the magnetic domain (default: ``None``).
    boundary_nodes : int, optional
        Number of boundary nodes for boundary value computations (default: ``1``). 
        Can be 1, 2, or 3. In setups with multiple material domains, this value 
        cannot exceed the number of boundary cells in the same domain, counting 
        from the interface for each interface individually.
    iteration_depth : int, optional
        Iteration depth for the strain jump conditions (default: ``1``).
    """

    def __init__(self, 
                 terms,                     # magnetic field terms
                 solver = RKF45,            # used solver
                 no_precession = False,     # remove the precessional term from the LLG
                 mask_elastic = None,       # torch.tensor (bool), indicates where ud and pd are updated
                 C_sym = None,              # 2D list, 6 times 6, entries of 1 indicate that the corresponding stiffness tensor component is non-zero. 0 indicates that it is zero.
                 magnetic_x_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in x-direction
                 magnetic_y_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in y-direction
                 magnetic_z_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in z-direction
                 boundary_nodes = 1,        # number of boundary nodes used for boundary value computations
                 iteration_depth = 1,       # iteration depth for strain jump conditions
                 **kwargs):

        # field terms for llg time integration
        self._terms = terms
        self._no_precession = no_precession

        # ode solver
        self._solver = solver(self.dv, **kwargs)
         
        # mask for state.pd and state.ud update
        self._additional_mask_elastic = mask_elastic

        # limit of the magnetic domain
        if magnetic_x_limits == None:
            self._magnetic_x_limits = [0, None]
        else:
            self._magnetic_x_limits = magnetic_x_limits
        if magnetic_y_limits == None:
            self._magnetic_y_limits = [0, None]
        else:
            self._magnetic_y_limits = magnetic_y_limits
        if magnetic_z_limits == None:
            self._magnetic_z_limits = [0, None]
        else:
            self._magnetic_z_limits = magnetic_z_limits

        ix0, ix1 = self._magnetic_x_limits
        iy0, iy1 = self._magnetic_y_limits
        iz0, iz1 = self._magnetic_z_limits

        self.slice_m = tuple([slice(ix0, ix1), slice(iy0, iy1), slice(iz0, iz1)])

        # mask for the stiffness matrix
        if isinstance(C_sym, str):
            if C_sym == "cubic" or C_sym == "isotropic":
                self._C_mask = [[1, 1, 1, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0, 1]]
            else:
                raise Exception(self.__class__.__name__ + ": C_sym '" + C_sym + "' is unknown. C_sym can be 'isotropic' 'cubic' or a mask of shape (6,6)")
        elif C_sym is None:
            self._C_mask = []
            self._C_mask.append(6*[1])
            for i in range(5):
                self._C_mask.append(self._C_mask[-1])
        else:
            if isinstance(C_sym, torch.Tensor):
                C_sym = C_sym.detach().cpu().numpy()
            if isinstance(C_sym, (list, np.ndarray)):
                C_sym = np.array(C_sym)
            if not isinstance(C_sym, np.ndarray) or C_sym.shape != (6, 6):
                raise Exception(self.__class__.__name__+": C_sym was not successfully converted into a (6,6) array")
            self._C_mask = C_sym

        self.iteration_depth = iteration_depth

        if (boundary_nodes < 1) or (boundary_nodes > 3):
            raise Exception(self.__class__.__name__ + ": boundary_nodes must be 1, 2 or 3.")
        self._gradient_second_order_boundary = boundary_nodes > 2 
        self._bc_second_order = boundary_nodes > 1

        # neumann bcs
        self._neumann_bcs = []

        # obtain the individual terms that make up the stress matrix and the force components
        # the former is organized as a 3-dimensional list:
        #       the first two dimensions being the matrix components of the stress matrix
        #       the last dimension contains the individual terms making up the corresponding stress component (lambda functions)
        # the latter is organized as follows: 
        #       the first dimension contains the components of the force, i.e. 0 = x, 1 = y, 2 = z,
        #       such that terms that contribute to fx are found in self._f_terms[0]
        #       the second dimension separates the components by their second derivative, 
        #       i.e. 0 = terms with second derivative in x-direction (e.g. Dx[C[0,1]*Dy[u[1]]]), ...
        #       the last dimension contains the individual terms (lambda functions)

        self._sig_terms, self._f_terms, self._fm_terms = self._compile_terms()

    # -----------------------------------
    # Definition of the solution variable
    # -----------------------------------

    #@torch.compile  
    def _get_solution_variables(self, state):
        n = state.mesh.n
        v = torch.zeros((n+(9,)))

        v[:,:,:,:3] = state.m 
        v[:,:,:,6:] = state.pd

        ud = state.ud
        if hasattr(state, "bcs"):  
            if "ud" in state.bcs:
                for bc in state.bcs["ud"]:
                    ud[bc.mask] = torch.flatten(bc.condition(state), end_dim=-2)

        v[:,:,:,3:6] = ud

        #v.detach().cpu().numpy().reshape(-1, order = 'F')

        return v
    
    #@torch.compile  
    def _set_solution_variables(self, state, v):
        # v = state.Tensor(v.reshape(state.mesh.n + (9,), order = "F"))

        state.m = v[:,:,:,:3]
        state.ud = v[:,:,:,3:6]
        state.pd = v[:,:,:,6:]

        if hasattr(state, "bcs"):
            if "ud" in state.bcs:
                for bc in state.bcs["ud"]:
                    state.ud[bc.mask] = torch.flatten(bc.condition(state), end_dim=-2)
    
    # ------------------------------
    # Definition of the problem area
    # ------------------------------

    def _get_mask_elastic(self, state):
        mask = self._additional_mask_elastic
        n = state.mesh.n

        if mask == None:
            mask = torch.ones(n, dtype=torch.bool)

        if hasattr(state, "bcs"):
            if "ud" in state.bcs:
                for bc in state.bcs["ud"]:
                    mask = torch.logical_and(mask, torch.logical_not(bc.mask))
        
        elastic_mask = torch.zeros(n+(3,), dtype=torch.bool)
        elastic_mask[:,:,:,0] = mask
        elastic_mask[:,:,:,1] = mask
        elastic_mask[:,:,:,2] = mask

        return elastic_mask
    
    # ----------------------------------------------------
    # Construction of the PDE RHS and stress (sigma) terms
    # ----------------------------------------------------
    
    def _compile_terms(self):
        compiler = DerivTermCompiler()

        # -------------
        # Collect terms
        # -------------

        # Initialize all strain terms
        # ... in Voigt notation
        eps = []
        eps.append([EpsTerm(i_u=0, i_x=0)])
        eps.append([EpsTerm(1,1)])
        eps.append([EpsTerm(2,2)])
        eps.append([EpsTerm(1,2), EpsTerm(2,1)])
        eps.append([EpsTerm(0,2), EpsTerm(2,0)])
        eps.append([EpsTerm(0,1), EpsTerm(1,0)])

        eps_m = []
        eps_m.append(EpsMTerm(0,0))
        eps_m.append(EpsMTerm(1,1))
        eps_m.append(EpsMTerm(2,2))
        eps_m.append(EpsMTerm(1,2))
        eps_m.append(EpsMTerm(0,2))
        eps_m.append(EpsMTerm(0,1))

        # multiply them by the stiffness tensor and obtain the stress terms
        sig_v = []
        sig_m_v = []
        for i in range(6):
            terms = []
            terms_m = []
            for j in range(6):
                # if the user informed the solver about the symmetry of the problem, exclude terms with a zero-valued stiffness tensor component
                if self._C_mask[i][j] == 0:
                    logging.info_green("Warning: due to the selected symmetry of C, a term with C_{"+str(i)+","+str(j)+"} was excluded!")
                else:
                    terms_m.append(eps_m[j].multiply_Cij(i,j))
                    for eps_term in eps[j]:
                        terms.append(eps_term.multiply_Cij(i,j))

            sig_v.append(terms)
            sig_m_v.append(terms_m)

        # rearrange the stress terms from Voigt notation into matrix form
        sig = []
        sig.append([sig_v[0], sig_v[5], sig_v[4]])
        sig.append([sig_v[5], sig_v[1], sig_v[3]])
        sig.append([sig_v[4], sig_v[3], sig_v[2]])

        sig_m = []
        sig_m.append([sig_m_v[0], sig_m_v[5], sig_m_v[4]])
        sig_m.append([sig_m_v[5], sig_m_v[1], sig_m_v[3]])
        sig_m.append([sig_m_v[4], sig_m_v[3], sig_m_v[2]])

        # -------------
        # compile terms
        # -------------

        sig_lambdas = []
        f_lambdas = []
        fm_lambdas = []
        for i in range(3):
            sig_row = []
            f_row = []
            fm_row = []
            for j in range(3):
                sig_element = []
                f_element = []
                fm_element = []
                for sig_term in sig[i][j]:
                    term_function = compiler.compile(self, sig_term)
                    sig_element.append(term_function)

                    term = sig_term.differentiate(j)
                    term_function = compiler.compile(self, term)
                    f_element.append(term_function)

                for sig_m_term in sig_m[i][j]:
                    term = sig_m_term.differentiate(j)
                    term_function = compiler.compile(self, term)
                    fm_element.append(term_function)
                
                sig_row.append(sig_element)
                f_row.append(f_element)
                fm_row.append(fm_element)

            sig_lambdas.append(sig_row)
            f_lambdas.append(f_row)
            fm_lambdas.append(fm_row)

        return sig_lambdas, f_lambdas, fm_lambdas

    # ----------------------------------------------------
    # Definition of 1st, 2nd, and mixed derivatives
    # ----------------------------------------------------

    class DiffData:
        # Object that holds forward differences and first derivatives obtained from the midpoint rule
        def __init__(self):
            self._gradient_ud = []
            self._gradient_m = []

        def set_expanded_dx(self, state):
            dx_exp, dy_exp, dz_exp = state.mesh.dx_tensor

            dx_exp = dx_exp.unsqueeze(1).unsqueeze(2)
            dy_exp = dy_exp.unsqueeze(0).unsqueeze(2)
            dz_exp = dz_exp.unsqueeze(0).unsqueeze(0)

            dx_exp = dx_exp.expand_as(torch.zeros(state.mesh.n))
            dy_exp = dy_exp.expand_as(torch.zeros(state.mesh.n))
            dz_exp = dz_exp.expand_as(torch.zeros(state.mesh.n))
            self.dx_epx = dx_exp, dy_exp, dz_exp

        def set_gradient_ud(self, grad_ux, grad_uy, grad_uz):
            data = []
            data.append(grad_ux)
            data.append(grad_uy)
            data.append(grad_uz)

            self._gradient_ud = data

        def set_gradient_m(self, grad_mx, grad_my, grad_mz):
            data = []
            data.append(grad_mx)
            data.append(grad_my)
            data.append(grad_mz)

            self._gradient_m = data

        def set_B_jump_conditions(self, Bl, Br):
            self._Bl_jump_conditions = Bl
            self._Br_jump_conditions = Br

        @property
        def gradient_ud(self):
            return self._gradient_ud
        
        @property
        def gradient_m(self):
            return self._gradient_m
    
    def _2nd_derivative(self, state, i_u, i_x, ij_C):
        if (state.mesh.pbc[i_x] == 0) and (state.mesh.n[i_x]>1):
            return self._2nd_derivative_homogeneous_Neumann(state, i_u, i_x, ij_C)
        else:
            return self._2nd_derivative_with_pbc(state, i_u, i_x, ij_C)

    #@torch.compile
    def _2nd_derivative_homogeneous_Neumann(self, state, i_u, i_x, ij_C):
        a = torch.zeros(state.mesh.n)

        dx = self.diff_data.dx_epx[i_x]

        slice_0, slice_1 = _get_diff_slices(3, i_x)

        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        C_denom = (dx[slice_0]*C[slice_1] + dx[slice_1]*C[slice_0])
        C_avg = 2.*C[slice_0]*C[slice_1] / C_denom
        C_avg.nan_to_num(posinf=0, neginf=0) # C could be 0 if not proper C_mask is set

        # Look up: [Bx_xl, Bx_yl, Bx_zl], [By_xl, By_yl, By_zl], [Bz_xl, Bz_yl, Bz_zl]
        Bl = self.diff_data._Bl_jump_conditions[i_u][i_x]
        Br = self.diff_data._Br_jump_conditions[i_u][i_x]
        jump = -(Br[slice_0] - Bl[slice_1]) / C_denom
        jump.nan_to_num(posinf=0, neginf=0)

        diff = state.ud[slice_1+(i_u,)] - state.ud[slice_0+(i_u,)] 

        a[slice_0] += C_avg*diff # (ud_i+1 - ud_i) * (C_i+1*dx_i + dx_1+i*C_i) / (2*C_i+1*C_i)
        a[slice_1] -= C_avg*diff # -(ud_i - ud_i-1) * (dx_i-1*C_i + dx_i*C_i-1)

        a[slice_0] += dx[slice_1]*C[slice_0]*jump
        a[slice_1] += dx[slice_0]*C[slice_1]*jump

        set_slices = [slice(None)]*dx.ndim 
        set_slices[i_x] = slice(1, -1)
        set_slices = tuple(set_slices)

        l_slices = [slice(None)]*dx.ndim 
        l_slices[i_x] = slice(2,None)
        l_slices = tuple(l_slices)

        r_slices = [slice(None)]*dx.ndim 
        r_slices[i_x] = slice(None,-2)
        r_slices = tuple(r_slices)

        dx_denom = torch.clone(dx)
        dx_denom[set_slices] += 0.5 * (dx[l_slices]+dx[r_slices])
        dx_denom[set_slices] /= 2.
        
        return (a/(dx_denom)).nan_to_num(posinf=0, neginf=0)
    
    #@torch.compile
    def _2nd_derivative_with_pbc(self, state, i_u, i_x, ij_C):
        # second order accurate for regular and only first order accurate for irregular grids
        a = torch.zeros(state.mesh.n)

        dx = self.diff_data.dx_epx[i_x]
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        C_next = torch.roll(C, -1, dims=i_x) # positive shift, to align this with the definition of the forward differences
        dx_next = torch.roll(dx, -1, dims=i_x) # positive shift, to align this with the definition of the forward differences
        C_denom = (C_next*dx + C*dx_next) # at 0: C_1*dx_0 + C_0*dx_1, at N: C_0*dx_N + C_N*dx_0
        C_avg = 2.*C_next * C / C_denom
        C_avg.nan_to_num(posinf=0, neginf=0) # C could be 0 if not proper C_mask is set

        diff = torch.roll(state.ud[...,i_u],-1,i_x) - state.ud[...,i_u]
        a += C_avg * diff
        a -= torch.roll(a, +1, i_x)

        # Look up: [Bx_xl, Bx_yl, Bx_zl], [By_xl, By_yl, By_zl], [Bz_xl, Bz_yl, Bz_zl]
        Bl = self.diff_data._Bl_jump_conditions[i_u][i_x]
        Br = self.diff_data._Br_jump_conditions[i_u][i_x]

        jump = -(Br - torch.roll(Bl, -1, i_x)) / C_denom
        jump.nan_to_num(posinf=0, neginf=0) # C could be 0 if not proper C_mask is set
        
        a += torch.roll(dx, -1, i_x)*C*jump 
        a += torch.roll(dx, +1, i_x)*C*torch.roll(jump, +1, i_x)
        
        dx_denom = torch.clone(dx)
        dx_denom += 0.5*(torch.roll(dx, +1, i_x) + torch.roll(dx, -1, i_x))
        dx_denom /= 2.

        return (a/dx_denom).nan_to_num(posinf=0, neginf=0)
    
    #@torch.compile
    def _mixed_derivative(self, state, i_u, i_x1, i_x2, ij_C):
        diff_1st = self.diff_data.gradient_ud[i_u][i_x1]
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        # Note: using C jump conditions here allows to handle jumps to vaccum correctly
        a = gradient_with_pbc(C*diff_1st, state.mesh, dim=i_x2, C=[C], second_order_boundary=self._gradient_second_order_boundary)[0]

        return a
    
    #@torch.compile
    def _weighted_first_derivative(self, state, i_u, i_x, ij_C):
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]
        return C * self.diff_data.gradient_ud[i_u][i_x]
    
    #@torch.compile 
    def _main_diag_sigM_derivative(self, state, i_m, i_x, ij_C):
        slice_m = self.slice_m
        l100 = state.material["lambda_100"][slice_m+(0,)]
        C = state.material["C"][slice_m+(ij_C[0],ij_C[1])]
        mi = state.m[slice_m+(i_m,)]
        dmi = self.diff_data.gradient_m[i_m][i_x]

        return 3.*l100*C*mi*dmi 
    
    #@torch.compile
    def _off_diag_sigM_derivative(self, state, i_m, j_m, i_x, ij_C):
        slice_m = self.slice_m
        l111 = state.material["lambda_111"][slice_m+(0,)]
        C = state.material["C"][slice_m+(ij_C[0],ij_C[1])]
        mi = state.m[slice_m+(i_m,)]
        mj = state.m[slice_m+(j_m,)]
        dmi = self.diff_data.gradient_m[i_m][i_x]
        dmj = self.diff_data.gradient_m[j_m][i_x]

        return 3.*l111*C*(mi*dmj + mj*dmi)

    # ----------------------------------------------------
    # Helpers
    # ----------------------------------------------------

    def _update_diff_data(self, state):
        diff_data = self.DiffData()
        diff_data.set_expanded_dx(state)

        # Get 2nd order 1st derivatives of m
        A = state.material["A"][...,0]

        slc_m = self.slice_m
        grad_mx = gradient_with_pbc(state.m[...,0], state.mesh, [0,1,2], [A, A, A], slices=slc_m, second_order_boundary=self._gradient_second_order_boundary)
        grad_my = gradient_with_pbc(state.m[...,1], state.mesh, [0,1,2], [A, A, A], slices=slc_m, second_order_boundary=self._gradient_second_order_boundary)
        grad_mz = gradient_with_pbc(state.m[...,2], state.mesh, [0,1,2], [A, A, A], slices=slc_m, second_order_boundary=self._gradient_second_order_boundary)

        diff_data.set_gradient_m(grad_mx, grad_my, grad_mz)

        # Get 2nd order 1st derivatives of u 
        mxl = torch.clone(state.m)
        mxr = torch.clone(state.m)
        myl = torch.clone(state.m)
        myr = torch.clone(state.m)
        mzl = torch.clone(state.m)
        mzr = torch.clone(state.m)

        dx_exp = state.mesh.dx_tensor[0].reshape(-1,1,1,1)
        dy_exp = state.mesh.dx_tensor[1].reshape(1,-1,1,1)
        dz_exp = state.mesh.dx_tensor[2].reshape(1,1,-1,1)

        dmdx = torch.zeros(state.mesh.n+(3,))
        dmdx[slc_m+(0,)] = grad_mx[0]
        dmdx[slc_m+(1,)] = grad_my[0]
        dmdx[slc_m+(2,)] = grad_mz[0]

        dmdy = torch.zeros(state.mesh.n+(3,))
        dmdy[slc_m+(0,)] = grad_mx[1]
        dmdy[slc_m+(1,)] = grad_my[1]
        dmdy[slc_m+(2,)] = grad_mz[1]

        dmdz = torch.zeros(state.mesh.n+(3,))
        dmdz[slc_m+(0,)] = grad_mx[2]
        dmdz[slc_m+(1,)] = grad_my[2]
        dmdz[slc_m+(2,)] = grad_mz[2]

        mxl -= 0.5*dmdx*dx_exp
        mxr += 0.5*dmdx*dx_exp
        myl -= 0.5*dmdy*dy_exp
        myr += 0.5*dmdy*dy_exp
        mzl -= 0.5*dmdz*dz_exp
        mzr += 0.5*dmdz*dz_exp

        m_data = [mxl, mxr], [myl, myr], [mzl, mzr]
        C = _get_C_jump_conditions(state)
        Bl_sigM, Br_sigM = _get_sigM_jump_conditions(state, m_data)

        grad_ud_x = gradient_with_pbc(state.ud[...,0], state.mesh, [0,1,2], C[0], Bl_sigM[0], Br_sigM[0], second_order_boundary=self._gradient_second_order_boundary)
        grad_ud_y = gradient_with_pbc(state.ud[...,1], state.mesh, [0,1,2], C[1], Bl_sigM[1], Br_sigM[1], second_order_boundary=self._gradient_second_order_boundary)
        grad_ud_z = gradient_with_pbc(state.ud[...,2], state.mesh, [0,1,2], C[2], Bl_sigM[2], Br_sigM[2], second_order_boundary=self._gradient_second_order_boundary)

        diff_data.set_B_jump_conditions(Bl_sigM, Br_sigM)

        if (self.iteration_depth > 0):
            for iter in range(self.iteration_depth):
                gradient_data = grad_ud_x, grad_ud_y, grad_ud_z

                Bl_eps, Br_eps = _get_B_jump_conditions(state, gradient_data)

                Bxl = [Bl_eps[0][0] + Bl_sigM[0][0], Bl_eps[0][1] + Bl_sigM[0][1], Bl_eps[0][2] + Bl_sigM[0][2]]
                Byl = [Bl_eps[1][0] + Bl_sigM[1][0], Bl_eps[1][1] + Bl_sigM[1][1], Bl_eps[1][2] + Bl_sigM[1][2]]
                Bzl = [Bl_eps[2][0] + Bl_sigM[2][0], Bl_eps[2][1] + Bl_sigM[2][1], Bl_eps[2][2] + Bl_sigM[2][2]]

                Bxr = [Br_eps[0][0] + Br_sigM[0][0], Br_eps[0][1] + Br_sigM[0][1], Br_eps[0][2] + Br_sigM[0][2]]
                Byr = [Br_eps[1][0] + Br_sigM[1][0], Br_eps[1][1] + Br_sigM[1][1], Br_eps[1][2] + Br_sigM[1][2]]
                Bzr = [Br_eps[2][0] + Br_sigM[2][0], Br_eps[2][1] + Br_sigM[2][1], Br_eps[2][2] + Br_sigM[2][2]]

                Bl = Bxl, Byl, Bzl 
                Br = Bxr, Byr, Bzr 

                grad_ud_x[:] = gradient_with_pbc(state.ud[...,0], state.mesh, [0,1,2], C[0], Bl[0], Br[0], second_order_boundary=self._gradient_second_order_boundary)[:]
                grad_ud_y[:] = gradient_with_pbc(state.ud[...,1], state.mesh, [0,1,2], C[1], Bl[1], Br[1], second_order_boundary=self._gradient_second_order_boundary)[:]
                grad_ud_z[:] = gradient_with_pbc(state.ud[...,2], state.mesh, [0,1,2], C[2], Bl[2], Br[2], second_order_boundary=self._gradient_second_order_boundary)[:]
                
                diff_data.set_B_jump_conditions(Bl, Br)

        diff_data.set_gradient_ud(grad_ud_x, grad_ud_y, grad_ud_z)

        self.diff_data = diff_data

    # ----------------------------------------------------
    # PDE RHS setters
    # ----------------------------------------------------
    
    @timedmethod
    def get_fm(self, state):
        # Returns the magnetic part of the force field
        m = state.m
        slice_m = self.slice_m
        
        """ magnetic strain forces """
        mesh = state.mesh
        n = mesh.n
        f_m = torch.zeros((n+(3,3)))

        for i in range(3):
            for j in range(3):
                for term in self._fm_terms[i][j]:
                    f_m[slice_m+(i,j)] += term(state)
        return f_m
    
    def _get_boundary_f(self, t0, s1, s2, hl, hr): 
        # This is actually the midpoint, since the boundary value is outside the node grid
        f_bdr = s2*hl**2. - t0*hr**2. + (hr**2.-hl**2.)*s1
        f_bdr /= hr*hl**2. + hl*hr**2.
        return f_bdr

    #@timedmethod
    def dpd(self, state):
        """ constants """
        n = state.mesh.n

        """ get first derivatives """
        self._update_diff_data(state)

        """ get magnetic stress """
        eps_m = epsilon_m(state)
        sig_m = sigma(state, eps_m)

        """ get stress """
        sig_ii = torch.zeros(n + (3,))
        for i in range(3):
            for sig_term in self._sig_terms[i][i]:
                sig_ii[:,:,:,i] += sig_term(state)

        
        sig_ij = torch.zeros(n + (3,))
        for sig_term in self._sig_terms[1][2]:
            sig_ij[:,:,:,0] += sig_term(state)

        for sig_term in self._sig_terms[0][2]:
            sig_ij[:,:,:,1] += sig_term(state)

        for sig_term in self._sig_terms[0][1]:
            sig_ij[:,:,:,2] += sig_term(state)

        sig_ii -= sig_m[...,:3]
        sig_ij -= sig_m[...,3:]

        """ get bulk forces """
        f_ij = torch.zeros(n + (3,3))

        for i in range(3):
            for j in range(3):
                for term in self._f_terms[i][j]:
                    f_ij[:,:,:,i,j] += term(state)

        """ get expanded dx """
        dx1_exp =  state.mesh.dx_tensor[0].unsqueeze(1).unsqueeze(2).expand_as(state.ud[...,0])
        dx2_exp =  state.mesh.dx_tensor[1].unsqueeze(0).unsqueeze(2).expand_as(state.ud[...,0])
        dx3_exp =  state.mesh.dx_tensor[2].unsqueeze(0).unsqueeze(0).expand_as(state.ud[...,0])
        dx_exp = dx1_exp, dx2_exp, dx3_exp

        """ set boundary conditions """
        ft = torch.zeros(state.mesh.n + (3,3))
        ft_weigth = torch.zeros(state.mesh.n + (3,3), dtype=int)
        for t_bc in self._neumann_bcs:
            t_dim = t_bc.plane.dim
            t_trans_dim1 = t_bc.plane.trans_dim1
            t_trans_dim2 = t_bc.plane.trans_dim2
            t_sign = t_bc.plane.sign
            t_val = t_bc.condition(state)
            t_slc = t_bc.plane.slices

            if self._bc_second_order and n[t_dim] > 1:
                if self._gradient_second_order_boundary:
                    hl = 0.5*dx_exp[t_dim][t_slc] # distance to from s1 to boundary
                    hr = hl + 0.5*(torch.roll(dx_exp[t_dim], t_sign, dims=t_dim)[t_slc]) # distance between s1 and s2

                    # t bc: main diagonal
                    t0 = t_val[...,t_dim]
                    s1 = sig_ii[t_slc+(t_dim,)] # sig value on node
                    s2 = torch.roll(sig_ii[...,t_dim], t_sign, dims=t_dim)[t_slc] # sig value on next node            

                    g_i = self._get_boundary_f(t0, s1, s2, hl, hr)
                    g_i *= -t_sign
                    
                    # t bc: remaining out of plane derivatives
                    t0 = t_val[...,t_trans_dim1]
                    s1 = sig_ij[t_slc+(t_trans_dim2,)] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                    s2 = torch.roll(sig_ij[...,t_trans_dim2], t_sign, dims=t_dim)[t_slc] # Note: sig_ij is counted "inverse", as yz, xz, xy 

                    g_j = self._get_boundary_f(t0, s1, s2, hl, hr)
                    g_j *= -t_sign

                    t0 = t_val[...,t_trans_dim2]
                    s1 = sig_ij[t_slc+(t_trans_dim1,)] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                    s2 = torch.roll(sig_ij[...,t_trans_dim1], t_sign, dims=t_dim)[t_slc] # Note: sig_ij is counted "inverse", as yz, xz, xy 

                    g_k = self._get_boundary_f(t0, s1, s2, hl, hr)
                    g_k *= -t_sign

                else:
                    h = dx_exp[t_dim][t_slc] + 0.5*(torch.roll(dx_exp[t_dim], t_sign, dims=t_dim)[t_slc])

                    # t bc: main diagonal
                    s_bdr = torch.roll(sig_ii[...,t_dim], t_sign, dims=t_dim)[t_slc]
                    g_i = t_sign*(t_val[:,:,t_dim] - s_bdr)/h # compute value
                    
                    # t bc: remaining out of plane derivatives
                    s_bdr = torch.roll(sig_ij[...,t_trans_dim2], t_sign, dims=t_dim)[t_slc] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                    g_j = t_sign*(t_val[:,:,t_trans_dim1] - s_bdr)/h
                    
                    s_bdr = torch.roll(sig_ij[...,t_trans_dim1], t_sign, dims=t_dim)[t_slc] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                    g_k = t_sign*(t_val[:,:,t_trans_dim2] - s_bdr)/h

            else:
                # Note: When forward and backward difference are used in this scheme for the boundary values of sig,
                # they are better estimates for the value at the cell boundary (Mean Value Theorem). 
                # Thus, h is the full cell width dx here
                h = dx_exp[t_dim][t_slc]

                # t bc: main diagonal
                s_bdr = sig_ii[t_slc+(t_dim,)] # sig value on node
                g_i = t_sign*(t_val[:,:,t_dim] - s_bdr)/h # compute value
                
                # t bc: remaining out of plane derivatives
                s_bdr = sig_ij[t_slc+(t_trans_dim2,)] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                g_j = t_sign*(t_val[:,:,t_trans_dim1] - s_bdr)/h
                
                s_bdr = sig_ij[t_slc+(t_trans_dim1,)] # Note: sig_ij is counted "inverse", as yz, xz, xy 
                g_k = t_sign*(t_val[:,:,t_trans_dim2] - s_bdr)/h


            ft[t_slc+(t_dim, t_dim)] += g_i 
            ft[t_slc+(t_trans_dim1, t_dim)] += g_j 
            ft[t_slc+(t_trans_dim2, t_dim)] += g_k

            ft_weigth[t_slc+(t_dim, t_dim)] += 1
            ft_weigth[t_slc+(t_trans_dim1, t_dim)] += 1 
            ft_weigth[t_slc+(t_trans_dim2, t_dim)] += 1
        
        """ get forces due to magnetic strain """
        f_ij -= self.get_fm(state)

        """ set boundary conditions """
        bc_mask = ft_weigth > 0 
        f_ij[bc_mask] = ft[bc_mask] / ft_weigth[bc_mask]

        """ add up all force contributions """
        self.f_ij = f_ij
        f_el = f_ij.sum(dim=-1)

        eta = state.material["eta"]
        return (f_el - eta*state.pd) * self._get_mask_elastic(state)

    @timedmethod
    def dud(self, state):
        rho = state.material["rho"]

        return (state.pd / rho) * self._get_mask_elastic(state) 

    def dv(self, t, v, state, alpha = None):
        
        self._set_solution_variables(state, v)

        # LLG 
        dm = self.dm(t, state.m, state, alpha) # sets t and m before dm is calculated

        # LE
        dud = self.dud(state)
        dpd = self.dpd(state)

        # composite variable
        n = state.mesh.n
        dv = torch.zeros((n[0],n[1],n[2],9))

        dv[:,:,:,:3] = dm
        dv[:,:,:,3:6] = dud
        dv[:,:,:,6:] = dpd

        return dv

    def E(self, state):
        return sum([term.E(state) for term in self._terms])
    
    def U_el(self, state):
        eps = epsilon(state, second_order_boundary=self._gradient_second_order_boundary)
        eps_el = eps - epsilon_m(state)
        sig_el = sigma(state, eps_el)

        zeta = 0.5*eps_el*sig_el

        E_el =  zeta * state.mesh.cell_volumes
        return E_el.sum()
    
    def U(self, state):
        eps = epsilon(state, second_order_boundary=self._gradient_second_order_boundary)
        sig = sigma(state, eps)

        zeta = 0.5*eps*sig

        E_el =  zeta * state.mesh.cell_volumes
        return E_el.sum()

    def T_el(self, state):
        T = 0.5 * state.pd**2. / state.material["rho"]
        return (T*state.mesh.cell_volumes*self._get_mask_elastic(state)).sum()
    
    def _update_neumann_bcs(self, state):
        self._neumann_bcs = []

        def _t_add(state, t0, t1_cond, trans_slice1, trans_slice2):
            t0[trans_slice1, trans_slice2] += t1_cond(state)
            return t0 
        
        def _apply_t_additions(state, t0, addition_list):
            t_out = torch.clone(t0)
            for addition in addition_list:
                t_out = addition(state, t_out)
            return t_out
        
        t_additions = [[[],[]],[[],[]],[[],[]]] # Note: 3*[[[],[]]] due to call by reference not possible
        
        # check for Neumann bcs set by the user
        if hasattr(state, "bcs"):
            if "t" in state.bcs:
                for bc in state.bcs["t"]:
                    # check if there is a conflict with natural boundary conditions
                    pl = bc.plane
                    i_where = None
                    if pl.pos == 0:
                        assert pl.sign == -1
                        i_where = 0
                    elif pl.pos == -1 or pl.pos == state.mesh.n[pl.dim]-1:
                        assert pl.sign == 1
                        i_where = -1

                    # CASE: conflict, mark override
                    if i_where != None:
                        assert state.mesh.pbc[pl.dim] == 0
                        t_additions[pl.dim][i_where].append(lambda state, t0, bc=bc, pl=pl : _t_add(state, t0, bc.condition, pl.trans_slice1, pl.trans_slice2))
                    # CASE: no conflict, use bc as is
                    else:
                        self._neumann_bcs.append(bc)

        # set natural boundary conditions where no pbc are set
        for i in range(3):
            if state.mesh.pbc[i] == 0:
                plane_m = Plane(i, 0, -1)
                plane_p = Plane(i, -1, 1)

                # remove the dimension in which the plane lies from n
                n_bc = list(state.mesh.n)
                n_bc.pop(i)
                n_bc = tuple(n_bc)

                # set homogenous Neumann boundary conditions
                tm = torch.zeros(n_bc + (3,))
                tp = torch.zeros(n_bc + (3,))

                # check for overrides
                if (len(t_additions[i][0]) > 0):
                    cond_m = lambda state, t0=tm, additions=t_additions[i][0]: _apply_t_additions(state, t0, additions)
                else:
                    cond_m = lambda state, t0=tm : t0 

                if (len(t_additions[i][1]) > 0):
                    cond_p = lambda state, t0=tp, additions=t_additions[i][1] : _apply_t_additions(state, t0, additions)
                else:
                    cond_p = lambda state, t0=tp : t0

                # set boundaries
                bc_plane_m = PlaneBC(plane_m, cond_m)
                bc_plane_p = PlaneBC(plane_p, cond_p)
                
                self._neumann_bcs.append(bc_plane_m)
                self._neumann_bcs.append(bc_plane_p)

    @timedmethod
    def step(self, state, dt, rtol = 1e-5, atol = None, atol_m = 1e-5, atol_ud = 1e-15, atol_pd = 1e-2, **kwargs):
        self._update_neumann_bcs(state)
        v_in = self._get_solution_variables(state)

        if not isinstance(atol, torch.Tensor):
            atol = tuple(3*[atol_m] + 3*[atol_ud] + 3*[atol_pd])
            atol = torch.tensor(atol)[None,None,None,:]

        state.t, v_out = self._solver.step(state.t, v_in, dt, state=state, rtol=rtol, atol=atol, **kwargs)
        self._set_solution_variables(state, v_out)
        normalize(state.m)
        logging.info_blue("[LLG + LE] step: dt= %g  t=%g" % (dt, state.t))

    @timedmethod
    def relax(self, state, maxiter = 500, dm_tol = 1e2, dud_tol = 1e-4, dt = 1e-11, rtol=1e-5, atol_m = 1e-5, atol_ud = 1e-15, atol_pd = 1e-2):
        t0 = state.t

        eta_in = state.material["eta"]
        state.material["eta"] = state.Constant(1e10)

        atol = tuple(3*[atol_m] + 3*[atol_ud] + 3*[atol_pd])
        atol = torch.tensor(atol)[None,None,None,:]

        self._update_neumann_bcs(state)

        for i in range(maxiter):
            # step
            v_in = self._get_solution_variables(state)
            state.t, v_out = self._solver.step(state.t, v_in, dt, state=state, rtol=rtol, atol=atol, alpha = 1.0)
            self._set_solution_variables(state, v_out)

            # determine rate of change of m
            dm = self.dm(state.t, state.m, state=state, alpha = 1.0)
            dm = dm[self.slice_m + (slice(None,3),)].abs().max() / constants.gamma # use same scaling as within minimizer

            # determine rate of change of ud
            dud = state.pd / state.material["rho"]
            dud = dud[self._get_mask_elastic(state)].nan_to_num(posinf=0, neginf=0).abs().max()

            logging.info_blue("[LLG+LE] relax: i=%d t=%g |dm|=%g, |dud|=%g" % (i, state.t-t0, dm, dud))
            if dm < dm_tol and dud < dud_tol:
                logging.info_green("[LLG+LE] relax: Successfully converged (iter=%d, dm_tol = %g, dud_tol = %g)" % (i, dm_tol, dud_tol))
                state.t = t0
                state.material["eta"] = eta_in
                state.pd = state.Constant((0.,0.,0.))
                return True

        logging.warning("[LLG+LE] relax: Terminated after maxiter = %d (dm = %g, dm_tol = %g, dud = %g, dud_tol = %g)" % (maxiter, dm, dm_tol,  dud, dud_tol))
        state.t = t0
        state.material["eta"] = eta_in
        state.pd = state.Constant((0.,0.,0.))
        return False
