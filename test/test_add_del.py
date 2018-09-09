# test add_del for `Network-processes`
#
# Copyright (C) 2017 Peter Mann
# 
# This file is part of `Network_processes`, for epidemic network 
# analytical results using Python.
#
# `Network_processes` is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# `Network_processes` is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with `Network_processes`. If not, see <http://www.gnu.org/licenses/gpl.html>.

from network_processes import *
import unittest

class addition_deletionTest(unittest.TestCase):
    '''Tests for `addition_deletion` class in `add_del.py` using an
    epyc lab simulation environment.'''
    
    def setUp( self ):
        '''Set up the parameters.'''
        # set lab test parameters
        self._lab = epyc.Lab()
        
        # initialise the experimental parameters
        self._lab['time'] = range(0,10)  # integration time
        self._lab['N'] = 5000              # network size
        self._lab['k_max'] = 30            # maximum degree
        self._lab['kmean'] = 10            # mean degree
        self._lab['class_dimension'] = 1   # number of equations per class
        
        self._lab['poisson'] = False       # flag for Poisson distribution
        self._lab['delta'] = True          # flag for delta function
        self._repetitions = 1              # repetitions at each point in the parameter space
    
    def testDegreeDist( self ):
        ''''Test an non-zero degree distribution for 
        steady-state network dynamics.'''
        # instance class
        e = addition_deletion()
        # perform the experiments
        self._lab.runExperiment(epyc.RepeatedExperiment(e, self._repetitions))
        # extract the results
        rc = (self._lab.results())[0]
        
        # perform tests
        self.assertTrue(sum(rc[epyc.Experiment.RESULTS]['sol']) > 0)

