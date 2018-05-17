# Initialisation for `Network-processes`
#
# Copyright (C) 2017 Peter Mann
# 
# This file is part of `Network-processes`, for epidemic network 
# analytical results using Python.
#
# `Network-processes` is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# `Network-processes` is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with `Network-processes`. If not, see <http://www.gnu.org/licenses/gpl.html>.

'''

`Network-processes` is a small library that implements analytical techniques to study 
dynamical processes on networks. Epidemic spreading processes on complex networks is 
an important area of current research. The dynamical properties of a disease can be 
captured using heterogeneous mean field theory (HMF), while the static, final-state 
properties are best portrayed through the use of generating functions (GFs). This 
library enables the study of an SIR process using both HMF and GFs.


'''
from .hmf import HMF 
from .gfs import GFs
from .sto import STO
from .gillespie import Gillespie