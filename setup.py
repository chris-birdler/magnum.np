#!/usr/bin/env python
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

from distutils.core import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='magnumnp',
      version='1.0.9',
      description='magnum.np finite-difference package for the solution of micromagnetic problems',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Florian Bruckner',
      author_email='florian.bruckner@univie.ac.at',
      url='http://gitlab.com/magnum.np/magnum.np',
      project_urls = {'Documentation': 'https://magnum.np.gitlab.io/magnum.np/',
                      'Changelog': 'https://gitlab.com/magnum.np/magnum.np/blob/main/CHANGELOG'
      },
      packages=['magnumnp', 'magnumnp.common', 'magnumnp.field_terms', 'magnumnp.loggers', 'magnumnp.solvers', 'magnumnp.solvers.ode_solvers', 'magnumnp.utils'],
      install_requires = [
            'torch',
            'numpy',
            'scipy',
            'setproctitle',
            'pyvista',
            'xitorch',
            'torchdiffeq'
            ]
     )
