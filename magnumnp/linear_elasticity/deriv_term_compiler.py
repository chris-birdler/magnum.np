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

import torch 
import torch._dynamo 

__all__ = ["DerivTermCompiler",
           "EpsTerm",
           "SigTerm",
           "EpsMTerm",
           "SigMTerm",
           "ForceMComponentTerm",
           "ForceComponentTerm"]

# TODO: replace torch.compile by torch.jit.trace(f, example_tensor)
torch._dynamo.config.cache_size_limit = 1024
torch._dynamo.config.suppress_errors = True

class DerivTermCompiler:
    # this object maps stress and force terms on lambda functions that return the given term as torch.tensor

    def compile(self, llg_with_le_solver, term):
        if isinstance(term, ForceComponentTerm):
            return torch.compile(self._compile_f_term(llg_with_le_solver, term))
        elif isinstance(term, ForceMComponentTerm):
            return torch.compile(self._compile_fm_term(llg_with_le_solver, term))
        elif isinstance(term, SigTerm):
            return torch.compile(self._compile_sig_term(llg_with_le_solver, term))
        else:
            raise SystemError("The derivative compiler does not recognize the type of the term provided to it. The type of the term is " + type(term) + " but should be ForceComponentTerm, ForceMComponentTerm or SigTerm!")

    def _compile_sig_term(self, solver, term):
        return lambda state, _solver=solver, iu=term.i_u, ix=term.i_x, ij=term.ij_C : \
            _solver._weighted_first_derivative(state, iu, ix, ij)
    
    def _compile_f_term(self, solver, term):
        if term.i_x1 == term.i_x2:
            return lambda state, _solver=solver, iu=term.i_u, ix=term.i_x1, ij=term.ij_C : \
                _solver._2nd_derivative(state, iu, ix, ij)     
        else:
            return lambda state, _solver=solver, iu=term.i_u, ix1=term.i_x1, ix2=term.i_x2, ij=term.ij_C : \
                _solver._mixed_derivative(state, iu, ix1, ix2, ij)
        
    def _compile_fm_term(self, solver, term):
        if term.i_m == term.j_m:
            return lambda state, _solver=solver, im=term.i_m, ix=term.i_x, ij=term.ij_C : \
                _solver._main_diag_sigM_derivative(state, im, ix, ij)
        else:
            return lambda state, _solver=solver, im=term.i_m, jm=term.j_m, ix=term.i_x, ij=term.ij_C : \
                _solver._off_diag_sigM_derivative(state, im, jm, ix, ij)

class EpsTerm:
    # information on the strain tensor terms
    def __init__(self, i_u, i_x):
        self._iu = i_u
        self._ix = i_x

    # multiplication with a component of the stiffness tensor yields a term of the stress tensor
    def multiply_Cij(self,i,j):
        return SigTerm(self._iu, self._ix, [i,j])

    # index of the displacement component that is differentiated
    @property
    def i_u(self):
        return self._iu
    
    # corresponding dimension in which the displacement is differentiated
    @property
    def i_x(self):
        return self._ix

class SigTerm:
    # information on the stress tensor terms
    def __init__(self, i_u, i_x, ij_C):
        self._iu = i_u
        self._ix = i_x
        self._ijC = ij_C

    # differentiation of stress terms yields terms of the force field
    def differentiate(self, j_x):
        return ForceComponentTerm(self._iu, self._ix, j_x, self._ijC)
    
    # indices (list with two entries, i and j) of the stiffness matrix that weights the derivative
    @property
    def ij_C(self):
        return self._ijC
    
    # index of the displacement component that is differentiated
    @property
    def i_u(self):
        return self._iu
    
    # corresponding dimension in which the displacement is differentiated
    @property
    def i_x(self):
        return self._ix
    
class EpsMTerm:
    def __init__(self, i_m, j_m):
        self._i_m = i_m 
        self._j_m = j_m 
    
    def multiply_Cij(self,i,j):
        return SigMTerm(self._i_m, self._j_m, [i,j])

    @property
    def i_m(self):
        return self._i_m
    
    @property
    def j_m(self):
        return self._j_m

class SigMTerm:
    def __init__(self, i_m, j_m, ij_C):
        self._i_m = i_m 
        self._j_m = j_m 
        self._ij_C = ij_C

    def differentiate(self, j_x):
        return ForceMComponentTerm(self._i_m, self._j_m, j_x, self._ij_C)

    @property
    def i_m(self):
        return self._i_m
    
    @property
    def j_m(self):
        return self._j_m
    
    @property
    def ij_C(self):
        return self._ij_C
    
class ForceMComponentTerm:
    def __init__(self, i_m, j_m, i_x, ij_C):
        self._i_m = i_m
        self._j_m = j_m
        self._i_x = i_x
        self._ij_C = ij_C

    @property
    def i_m(self):
        return self._i_m
    
    @property
    def j_m(self):
        return self._j_m
    
    # first derivative
    @property
    def i_x(self):
        return self._i_x
    
    @property
    def ij_C(self):
        return self._ij_C

class ForceComponentTerm:
    # information on the force field terms
    def __init__(self, i_u, i_x1, i_x2, ij_C):
        self._iu = i_u
        self._ix1 = i_x1
        self._ix2 = i_x2
        self._ijC = ij_C

    # as in SigTerm
    @property
    def i_u(self):
        return self._iu
    
    # first derivative
    @property
    def i_x1(self):
        return self._ix1
    
    # second derivative
    @property
    def i_x2(self):
        return self._ix2
    
    # as in SigTerm
    @property
    def ij_C(self):
        return self._ijC
