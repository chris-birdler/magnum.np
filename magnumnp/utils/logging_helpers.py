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

__all__ = ["log_cumsum", "log_diff"]


def log_cumsum(**kwargs):
    ''' Cummulate values '''
    items = list(kwargs.items())
    if len(items) > 1:
        raise ValueError("only a single parameter allowed") 
    key, value = items[0]
    try:  
        globals()[key] += value
    except: 
        globals()[key] = value 
    return globals()[key]
    

def log_diff(**kwargs):
    ''' Discrete difference '''
    items = list(kwargs.items())
    if len(items) > 1:
        raise ValueError("only a single parameter allowed") 
    key, value = items[0]
    try:  
        x0 = globals()[key]
    except: 
        x0 = 0
    diff = value - x0
    globals()[key] = value
    return diff
    

