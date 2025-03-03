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

from magnumnp.common import logging, timedmethod
from . import LLGSolver
from magnumnp.linear_elasticity.bcs import Plane, PlaneBC
from magnumnp.linear_elasticity.deriv_term_compiler import *
from magnumnp.linear_elasticity.strain import epsilon, epsilon_el, epsilon_m
from magnumnp.linear_elasticity.stress import sigma
from magnumnp.linear_elasticity.utils import gradient_with_pbc, _get_diff_slices

from .ode_solvers import RKF45
import torch

__all__ = ["LLGWithLESolver"]

class LLGWithLESolver(LLGSolver):
    def __init__(self, 
                 terms,                     # magnetic field terms
                 solver = RKF45,            # used solver
                 no_precession = False,
                 mask_elastic = None,       # torch.tensor (bool), indicates where ud and pd are updated
                 C_sym = None,              # 2D list, 6 times 6, entries of 1 indicate that the corresponding stiffness tensor component is non-zero. 0 indicates that it is zero.
                 magnetic_x_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in x-direction
                 magnetic_y_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in y-direction
                 magnetic_z_limits = None,  # list of len 2, upper and lower limit of the magnetic domain in z-direction
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
        if C_sym is None:
            self._C_mask = []
            self._C_mask.append(6*[1])
            for i in range(5):
                self._C_mask.append(self._C_mask[-1])
        elif C_sym == "cubic" or C_sym == "isotropic":
            self._C_mask = [[1, 1, 1, 0, 0, 0],
                            [1, 1, 1, 0, 0, 0],
                            [1, 1, 1, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0],
                            [0, 0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0, 1]]
        else:
            self._C_mask = C_sym

        # neumann bcs
        self._neumann_bcs = []

        # obtain the individual terms that make up the stress matrix and the force components
        # the former is orgainized as a 3 dimensional list:
        #       the first two dimensions being the matrix component of the stress matrix
        #       the last dimension contains the indiviudal terms making up the corresponding stress component (lambda functions)
        # the latter is organized as follows: 
        #       the first dimension are the components of the force i.e. 0 = x, 1 = y, 2 = z such that terms that contribute to fx are found in self._f_terms[0]
        #       the second dimension seperates the components by their second derivative i.e. 0 = terms with second derivative in x-direction (e.g. Dx[C[0,1]*Dy[u[1]]]), ...
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
    
    def _get_mask_magnetic(self, state):
        if self._mask_magnetic == None:
            return torch.ones(state.mesh.n)
        else:
            return self._mask_magnetic
    

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
        # collect terms
        # -------------

        # initialize all strain terms
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
                if self._C_mask[i][j] != 1:
                    logging.info_green("Warning: due to the selected symmetry of C, a term with C_{"+str(i)+","+str(j)+"} was excluded!")
                else:
                    terms_m.append(eps_m[j].multiply_Cij(i,j))
                    for eps_term in eps[j]:
                        terms.append(eps_term.multiply_Cij(i,j))

            sig_v.append(terms)
            sig_m_v.append(terms_m)

        # rearange the stress terms from Voight notation into matrix form
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

                    term = sig_term.deriviate(j)
                    term_function = compiler.compile(self, term)
                    f_element.append(term_function)

                for sig_m_term in sig_m[i][j]:
                    term = sig_m_term.deriviate(j)
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
    # Definition of 1st, 2nd and mixed derivatives
    # ----------------------------------------------------

    class DiffData:
        # object that holds forward differences and first derivatives obtained from the midpoint rule
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
        slice_0 = [slice(None)] * 3
        slice_1 = [slice(None)] * 3
        slice_0[i_x] = slice(0, -1)
        slice_1[i_x] = slice(1, None)
        slice_0 = tuple(slice_0)
        slice_1 = tuple(slice_1)

        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        C_denom = (dx[slice_0]*C[slice_1] + dx[slice_1]*C[slice_0])
        C_avg = 2.*C[slice_0]*C[slice_1] / C_denom
        C_avg = torch.nan_to_num(C_avg) # C could be 0 if not proper C_mask is set

        slc_m = self.slice_m
        m1 = torch.clone(state.m)
        m2 = torch.clone(state.m)

        m1[slc_m+(0,)] -= 0.5*self.diff_data.gradient_m[0][i_x]*dx[slc_m]
        m1[slc_m+(1,)] -= 0.5*self.diff_data.gradient_m[1][i_x]*dx[slc_m]
        m1[slc_m+(2,)] -= 0.5*self.diff_data.gradient_m[2][i_x]*dx[slc_m]

        m2[slc_m+(0,)] += 0.5*self.diff_data.gradient_m[0][i_x]*dx[slc_m]
        m2[slc_m+(1,)] += 0.5*self.diff_data.gradient_m[1][i_x]*dx[slc_m]
        m2[slc_m+(2,)] += 0.5*self.diff_data.gradient_m[2][i_x]*dx[slc_m]

        eps_m1 = epsilon_m(state, m1)[...,ij_C[1]]
        eps_m2 = epsilon_m(state, m2)[...,ij_C[1]]

        jump = (C[slice_0]*eps_m2[slice_0] - C[slice_1]*eps_m1[slice_1]) / C_denom
        jump = torch.nan_to_num(jump) # C could be 0 if not proper C_mask is set

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
        
        return a/dx_denom

    #@torch.compile
    def _2nd_derivative_with_pbc(self, state, i_u, i_x, ij_C):
        # second order accurate for regular and only first order accurate for irregular grids
        a = torch.zeros(state.mesh.n)

        dx = self.diff_data.dx_epx[i_x]
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        C_next = torch.roll(C, -1, dims=i_x) # positive shift, to align with the definition of the forward differences
        dx_next = torch.roll(dx, -1, dims=i_x) # positive shift, to align with the definition of the forward differences
        C_denom = (C_next*dx + C*dx_next) # at 0: C_1*dx_0 + C_0*dx_1, at N: C_0*dx_N + C_N*dx_0
        C_avg = 2.*C_next * C / C_denom
        C_avg = torch.nan_to_num(C_avg) # C could be 0 if not proper C_mask is set

        diff = torch.roll(state.ud[...,i_u],-1,i_x) - state.ud[...,i_u]
        a += C_avg * diff
        a -= torch.roll(a, +1, i_x)

        slc_m = self.slice_m
        m1 = torch.clone(state.m)
        m2 = torch.clone(state.m)

        m1[slc_m+(0,)] -= 0.5*self.diff_data.gradient_m[0][i_x]*dx[slc_m]
        m1[slc_m+(1,)] -= 0.5*self.diff_data.gradient_m[1][i_x]*dx[slc_m]
        m1[slc_m+(2,)] -= 0.5*self.diff_data.gradient_m[2][i_x]*dx[slc_m]

        m2[slc_m+(0,)] += 0.5*self.diff_data.gradient_m[0][i_x]*dx[slc_m]
        m2[slc_m+(1,)] += 0.5*self.diff_data.gradient_m[1][i_x]*dx[slc_m]
        m2[slc_m+(2,)] += 0.5*self.diff_data.gradient_m[2][i_x]*dx[slc_m]

        s1 = -epsilon_m(state, m1)[...,ij_C[1]]
        s2 = -epsilon_m(state, m2)[...,ij_C[1]]

        jump = (torch.roll(C, -1, i_x)*torch.roll(s1, -1, i_x) - C*s2) / C_denom
        jump = torch.nan_to_num(jump) # C could be 0 if not proper C_mask is set

        a += torch.roll(dx, -1, i_x)*C*jump 
        a += torch.roll(dx, +1, i_x)*C*torch.roll(jump, +1, i_x)

        dx_denom = torch.clone(dx)
        dx_denom += 0.5*(torch.roll(dx, +1, i_x) + torch.roll(dx, -1, i_x))
        dx_denom /= 2.

        return a/dx_denom
    
    #@torch.compile
    def _mixed_derivative(self, state, i_u, i_x1, i_x2, ij_C):
        diff_1st = self.diff_data.gradient_ud[i_u][i_x1]
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]

        a = gradient_with_pbc(C*diff_1st, state.mesh, dim=i_x2)[0]

        return a
    
    #@torch.compile
    def _weighted_first_derivative(self, state, i_u, i_x, ij_C):
        C = state.material["C"][:,:,:,ij_C[0],ij_C[1]]
        return C * self.diff_data.gradient_ud[i_u][i_x]
    
    #@torch.compile 
    def _main_diag_sig_derivative(self, state, i_m, i_x, ij_C):
        slice_m = self.slice_m
        l100 = state.material["lambda_100"][slice_m+(0,)]
        C = state.material["C"][slice_m+(ij_C[0],ij_C[1])]
        mi = state.m[slice_m+(i_m,)]
        dmi = self.diff_data.gradient_m[i_m][i_x]

        return 3.*l100*C*mi*dmi 
    
    #@torch.compile
    def _off_diag_sig_derivative(self, state, i_m, j_m, i_x, ij_C):
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

        """ get 2nd order 1st derivatives of m """
        A = state.material["A"][...,0]

        slc_m = self.slice_m
        grad_mx = gradient_with_pbc(state.m[...,0], state.mesh, [0,1,2], [A, A, A], slices=slc_m)
        grad_my = gradient_with_pbc(state.m[...,1], state.mesh, [0,1,2], [A, A, A], slices=slc_m)
        grad_mz = gradient_with_pbc(state.m[...,2], state.mesh, [0,1,2], [A, A, A], slices=slc_m)

        diff_data.set_gradient_m(grad_mx, grad_my, grad_mz)

        """ get 2nd order 1st derivatives of u """
        C = state.material["C"]
        C11 = C[...,0,0]
        C22 = C[...,1,1]
        C33 = C[...,2,2]
        C44 = C[...,3,3]
        C55 = C[...,4,4]
        C66 = C[...,5,5]

        C15 = C[...,0,4]
        C16 = C[...,0,5]
        C24 = C[...,1,3]
        C26 = C[...,1,5]
        C35 = C[...,2,4]
        C34 = C[...,2,3]
        C45 = C[...,3,4]
        C46 = C[...,3,5]
        C56 = C[...,4,5]

        mxl = torch.clone(state.m)
        mxr = torch.clone(state.m)
        myl = torch.clone(state.m)
        myr = torch.clone(state.m)
        mzl = torch.clone(state.m)
        mzr = torch.clone(state.m)

        #dx1_exp =  state.mesh.dx_tensor[0].unsqueeze(1).unsqueeze(2).expand_as(state.m)
        #dx2_exp =  state.mesh.dx_tensor[1].unsqueeze(0).unsqueeze(2).expand_as(state.m)
        #dx3_exp =  state.mesh.dx_tensor[2].unsqueeze(0).unsqueeze(0).expand_as(state.m)

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

        # for x derivatives
        Cxx = C11+C15+C16
        Cyx = C16+C66+C56
        Czx = C15+C56+C55

        eps_m = epsilon_m(state, mxl)
        epsMXX = eps_m[...,0]
        epsMXY = eps_m[...,5]
        epsMXZ = eps_m[...,4]
        
        Bxx_xl = C11*epsMXX + C15*epsMXZ + C16*epsMXY
        Byx_xl = C16*epsMXX + C66*epsMXY + C56*epsMXZ
        Bzx_xl = C15*epsMXX + C56*epsMXY + C55*epsMXZ

        eps_m = epsilon_m(state, mxr)
        epsMXX = eps_m[...,0]
        epsMXY = eps_m[...,5]
        epsMXZ = eps_m[...,4]

        Bxx_xr = C11*epsMXX + C15*epsMXZ + C16*epsMXY
        Byx_xr = C16*epsMXX + C66*epsMXY + C56*epsMXZ
        Bzx_xr = C15*epsMXX + C56*epsMXY + C55*epsMXZ

        # for y derivatives
        Cxy = C66+C26+C46
        Cyy = C26+C22+C24
        Czy = C46+C24+C44

        eps_m = epsilon_m(state, myl)
        epsMYY = eps_m[...,1]
        epsMXY = eps_m[...,5]
        epsMYZ = eps_m[...,3] 

        Bxy_yl = C66*epsMXY + C26*epsMYY + C46*epsMYZ
        Byy_yl = C26*epsMXY + C22*epsMYY + C24*epsMYZ
        Bzy_yl = C46*epsMXY + C24*epsMYY + C44*epsMYZ

        eps_m = epsilon_m(state, myr)
        epsMYY = eps_m[...,1]
        epsMXY = eps_m[...,5]
        epsMYZ = eps_m[...,3] 

        Bxy_yr = C66*epsMXY + C26*epsMYY + C46*epsMYZ
        Byy_yr = C26*epsMXY + C22*epsMYY + C24*epsMYZ
        Bzy_yr = C46*epsMXY + C24*epsMYY + C44*epsMYZ

        # for z derivatives
        Cxz = C55+C45+C35
        Cyz = C45+C44+C34
        Czz = C35+C34+C33

        eps_m = epsilon_m(state, mzl)
        epsMZZ = eps_m[...,2]
        epsMXZ = eps_m[...,4]
        epsMYZ = eps_m[...,3] 

        Bxz_zl = C55*epsMXZ + C45*epsMYZ + C35*epsMZZ
        Byz_zl = C45*epsMXZ + C44*epsMYZ + C34*epsMZZ
        Bzz_zl = C35*epsMXZ + C34*epsMYZ + C33*epsMZZ

        eps_m = epsilon_m(state, mzr)
        epsMZZ = eps_m[...,2]
        epsMXZ = eps_m[...,4]
        epsMYZ = eps_m[...,3] 

        Bxz_zr = C55*epsMXZ + C45*epsMYZ + C35*epsMZZ
        Byz_zr = C45*epsMXZ + C44*epsMYZ + C34*epsMZZ
        Bzz_zr = C35*epsMXZ + C34*epsMYZ + C33*epsMZZ

        grad_ux = gradient_with_pbc(state.ud[...,0], state.mesh, [0,1,2], [Cxx, Cxy, Cxz], [-Bxx_xl, -Bxy_yl, -Bxz_zl], [-Bxx_xr, -Bxy_yr, -Bxz_zr])
        grad_uy = gradient_with_pbc(state.ud[...,1], state.mesh, [0,1,2], [Cyx, Cyy, Cyz], [-Byx_xl, -Byy_yl, -Byz_zl], [-Byx_xr, -Byy_yr, -Byz_zr])
        grad_uz = gradient_with_pbc(state.ud[...,2], state.mesh, [0,1,2], [Czx, Czy, Czz], [-Bzx_xl, -Bzy_yl, -Bzz_zl], [-Bzx_xr, -Bzy_yr, -Bzz_zr])

        diff_data.set_gradient_ud(grad_ux, grad_uy, grad_uz)

        self.diff_data = diff_data

    # ----------------------------------------------------
    # PDE RHS setters
    # ----------------------------------------------------
    
    @timedmethod
    def get_fm(self, state):
        # returns the magnetic part of the force field
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
        f_bdr = (s2-s1)*hl**2. + (s1-t0)*hr**2
        f_bdr /= hl*hr*(hl+hr)
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

        """ get bulk forces """
        f_ij = torch.zeros(n + (3,3))
        for i in range(3):
            for j in range(3):
                for term in self._f_terms[i][j]:
                    f_ij[:,:,:,i,j] += term(state)

        """ get forces and set boundary conditions """
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
        
        dx1_exp =  state.mesh.dx_tensor[0].unsqueeze(1).unsqueeze(2).expand_as(state.ud[...,0])
        dx2_exp =  state.mesh.dx_tensor[1].unsqueeze(0).unsqueeze(2).expand_as(state.ud[...,0])
        dx3_exp =  state.mesh.dx_tensor[2].unsqueeze(0).unsqueeze(0).expand_as(state.ud[...,0])
        dx_exp = dx1_exp, dx2_exp, dx3_exp

        ft = torch.zeros(state.mesh.n + (3,3))
        ft_weigth = torch.zeros(state.mesh.n + (3,3), dtype=int)
        for t_bc in self._neumann_bcs:
            t_dim = t_bc.plane.dim
            t_trans_dim1 = t_bc.plane.trans_dim1
            t_trans_dim2 = t_bc.plane.trans_dim2
            t_sign = t_bc.plane.sign
            t_val = t_bc.condition(state)
            t_slc = t_bc.plane.slices

            if n[t_dim] > 1:
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
                h = 0.5*dx_exp[t_dim][t_slc] # get distance to the boundary

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
        f_el = f_ij.sum(dim=-1)

        eta = state.material["eta"]
        return (f_el - eta*state.pd) * self._get_mask_elastic(state)

    @timedmethod
    def dud(self, state):
        rho = state.material["rho"]

        return state.pd * self._get_mask_elastic(state) / rho

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
        eps_el = epsilon_el(state)
        sig_el = sigma(state, eps_el)

        zeta = 0.5*eps_el*sig_el

        E_el =  zeta * state.mesh.cell_volumes
        return E_el.sum()
    
    def U(self, state):
        eps = epsilon(state)
        sig = sigma(state, eps)

        zeta = 0.5*eps*sig

        E_el =  zeta * state.mesh.cell_volumes
        return E_el.sum()

    def T_el(self, state):
        T = 0.5 * state.pd**2. / state.material["rho"]
        return (T*state.mesh.cell_volumes).sum()
    
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
                    # CASE: no conflict, use bc as found
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
    """

    def _update_neumann_bcs(self, state):
        self._neumann_bcs = []

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
                t0 = torch.zeros(n_bc + (3,))
                bc_plane_m = PlaneBC(plane_m, t0)
                bc_plane_m.mask = plane_m.get_mask(state) # TODO: This is very ugly!
                bc_plane_p = PlaneBC(plane_p, t0)
                bc_plane_p.mask = plane_p.get_mask(state) # TODO: This is very ugly!
                self._neumann_bcs.append(bc_plane_m)
                self._neumann_bcs.append(bc_plane_p)

        # check for Neumann bcs set by the user
        if hasattr(state, "bcs"):
            if "t" in state.bcs:
                self._neumann_bcs += state.bcs["t"]
                for bc in self._neumann_bcs:
                    bc.mask = bc.plane.get_mask(state) # TODO: This is very ugly!
    """

    @timedmethod
    def step(self, state, dt, rtol = 1e-5, atol = (1e-5, 1e-5, 1e-5, 1e-15, 1e-15, 1e-15, 1e-2, 1e-2, 1e-2), **kwargs):
        self._update_neumann_bcs(state)
        v_in = self._get_solution_variables(state)
        state.t, v_out = self._solver.step(state.t, v_in, dt, state=state, rtol=rtol, atol=atol, **kwargs)
        self._set_solution_variables(state, v_out)
        logging.info_blue("[LLG + LE] step: dt= %g  t=%g" % (dt, state.t))
