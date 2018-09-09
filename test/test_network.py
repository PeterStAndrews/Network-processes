# test network for `Network-processes`
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
import networkx

class sample_experiment0( NETWORK ):
    '''A sample experiment that subclasses `NETWORK` and tests its 
    methods return expected results.'''
    
    def __init__(self):
        super(sample_experiment0, self).__init__()
        
    def do( self, params ):
        '''Executes methods in `NETWORK` and returns them as results.
        :param params: experimental parameters'''
        # dict to store results 
        rc = dict()
        
        # grab a copy of the network
        g = self._network
        
        # compute the degree distribution
        Pk = self.degree_distribution(g)
        
        # compute the average degree
        ave_k = self.average_degree(Pk)
        
        rc['Pk'] = Pk
        rc['ave_k'] = ave_k
        rc['network'] = g
        
        return rc
        
        
class NetworkTest( unittest.TestCase ):
    '''Test class for `NETWORK` class in `network.py`.'''
    
    def testNETWORK( self ):
        '''Test methods in `NETWORK` return expected results.'''
        # set up epyc lab
        self._lab = epyc.Lab()
        
        # initialise parameters
        self._lab[NETWORK.N] = 5000
        self._lab[NETWORK.AVERAGE_K] = 5
        
        # repetitions at each point in the parameter space
        self._repetitions = 1
        
        # instance experiment class, run the experiment and extract results
        e = sample_experiment0()
        self._lab.runExperiment(epyc.RepeatedExperiment(e, self._repetitions))
        rc = (self._lab.results())[0]
        
        # assert non-zero degree distribution
        self.assertTrue(sum(rc[epyc.Experiment.RESULTS]['Pk'].values()) > 0)
        
        # assert computed average degree is within +-1.0 of parameter `AVERAGE_K`
        comp = rc[epyc.Experiment.RESULTS]['ave_k'] 
        self.assertTrue(4 <= comp and comp <= 6)
        
        # assertTrue that network order is less than or equal to `N`
        self.assertTrue(len(rc[epyc.Experiment.RESULTS]['network']) <= self._lab[NETWORK.N])
        
        # assertTrue that network has no degree-zero nodes 
        self.assertFalse(networkx.isolates(rc[epyc.Experiment.RESULTS]['network']))
        
