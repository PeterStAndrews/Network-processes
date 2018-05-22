
# Network base class
#
# Copyright (C) 2017 Peter Mann
# 
# This file is part of `NetworkProcesses`, for epidemic network 
# analytical results using Python.
#
# `NetworkProcesses` is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# `NetworkProcesses` is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with `NetworkProcesses`. If not, see <http://www.gnu.org/licenses/gpl.html>.

import networkx
import epyc


class NETWORK( epyc.Experiment ):
    '''The base class for generating a network using networkx and
    defining basic analysis tools for later use. Each method will 
    subclass this to model the desired process. '''
    
    N = 'N' # order of the network
    AVERAGE_K = 'kmean' # average degree s

    def __init__(self):
        super(NETWORK, self).__init__()
        
    def configure( self, params ):
        '''Create a "prototype" network and store it 
        for later use.
        :param params: the experimental parameters'''
        epyc.Experiment.configure(self, params)
        
        # create the prototype network
        N = params[self.N]
        kmean = params[self.AVERAGE_K] + 0.0
        g = networkx.erdos_renyi_graph(N, kmean / N)
        
        # remove degree-zero nodes
        ks = g.degree()
        k0s = [ i for i in ks.keys() if ks[i] == 0 ]
        g.remove_nodes_from(k0s)
        
        # remove self-loops
        g.remove_edges_from(g.selfloop_edges())
        
        # store it for later
        self._prototype = g
        
    def setUp( self, params ):
        '''Set up a working network for this run of the experiment.
        This is useful when performing lab experiments.
        :param params: the experimental parameters'''
        epyc.Experiment.setUp(self, params)
        self._network = self._prototype.copy()

    def tearDown( self ):
        '''Delete the current network.'''
        epyc.Experiment.tearDown(self)
        self._network = None

    def degree_distribution( self, g ):
        '''Computes the degree distribution of the network and stores
        as a dictionary {degree: P_k}.
        :param g: the network'''
        Pk = {}
        order = g.order()
        inv_order = 1./order
        for node in g.nodes_iter():
            k = g.degree(node)
            Pk[k] = Pk.get(k,0) + inv_order
        return Pk
    
    def average_degree( self, Pk ):
        '''Returns the average degree of the degree distribution Pk.
        :param Pk: degree distribution'''
        ave_k = 0
        for k in Pk.keys():
            ave_k += k*Pk[k]
        return ave_k
    
