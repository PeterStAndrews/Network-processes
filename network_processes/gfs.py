# Generating functions base class
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


import networkx
import epyc


class GFs( epyc.Experiment ):
    '''Base for conducting epidemic analytical generating functions on
    networks. This class initialises the network and defines common 
    methods that will be used frequently. Sub-classes will define the 
    details of the experimental sections and will overide `do` in all
    but the basic scenarios.
    
    This class is built over an epyc.Experiment base class. In the 
    `configure` function we create a network of size `N` and average
    degree `k_average`. The next section includes methods for the 
    epyc base and for the evaluation of network properties.
    The `do()` function is the actual experiment to be run, this section
    is customised to the needs of the experimenter.
    
    :References:
    -------------
    .. [1] Newman. (2002) 'The spread of epidemic disease on networks'
    .. [2] Miller. (2007) 'Epidemic size and probability in populations 
       with heterogeneous infectivity and susceptibility'
    .. [3] Funk & Jansen (2010) `Interacting epidemics on overlay networks`
    '''
    
    N = 'N' # order of the network
    AVERAGE_K = 'kmean' # average degree s
    T = 'T'
    
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
        
    def return_arg( self, *arg ):
        '''Returns the argument `x` of a generating function.  
        :param arg: the argument of the function '''
        return arg[0]
    
    def G_0_generating_function( self, Pk, ave_k, x, *args, **kwargs ):
        r'''Computes a `G_0` generating function as 
        
         .. math::
         
             \sum_{k=0}^\infty p_k(x)^k
        
        :param Pk: degree distribution
        :param ave_k: average degree
        :param x: argument function'''
        G_0 = 0
        for k in Pk.keys():
            G_0 += Pk[k]*x(*args, **kwargs)**k
        return G_0
       
    def G_1_generating_function( self, Pk, ave_k, x, *args, **kwargs ):
        r'''Computes the `G_1` generating function as
        
         .. math::
         
             \frac{1}{\langle k \rangle}\sum_{k=0}^\infty (k+1)p_{k+1}(x)^k
         
        :param Pk: degree distribution
        :param ave_k: average degree
        :param x: argument function'''
        summation = 0
        k_min = min(Pk, key=int)
        for k in Pk.keys():
            if k > k_min:
                summation += k*Pk[k]*x(*args, **kwargs)**(k-1)
        return (summation + 0.0)/ave_k
 
    def G_1_derivative_generating_function( self, Pk, ave_k, x, *args, **kwargs ):
        r'''Computes the `G'_1` generating function derivative as
        
         .. math::
         
             \frac{1}{\langle k \rangle}\sum_{k=0}^\infty k(k-1)(k+1)p_{k+1}(x)^{k-2}
        
        :param Pk: degree distribution
        :param ave_k: average degree
        :param x: argument function'''
        summation = 0
        k_min = min(Pk, key=int)
        for k in Pk.keys():
            if k > (k_min + 1):
                summation += (k-1)*k*Pk[k]*x(*args, **kwargs)**(k-2)
        return (summation + 0.0)/ave_k
    
    def super_critical_GC( self, T, Pk, ave_k, u ):
        r'''Computes the fraction (`S_1`) of the network occupied by the 
        giant component (`GC`) of disease 1 in the super-critical regime. 
        
        :param T: transmissibility 
        :param Pk: degree distribution
        :param ave_k: average degree
        :param u: probability that a neighbour fails to transmit ''' 
        x = 1 - T + T*u
        S_1 = 1 - self.G_0_generating_function(Pk, ave_k, self.return_arg, x)
        return S_1
    
    def make_z( self, T, ave_k, Pk ):
        r'''Computes `z` the probability that a neighbour infected a node as
        
        .. math::
            z = 1 - G_1(1 - zT)
        
        :param T: Transmissibility 
        :param ave_k: average degree
        :param Pk: degree distribution'''
        z = 1
        for i in range(1000):
            z = 1 - z*T
            z = 1 - self.G_1_generating_function(Pk, ave_k, self.return_arg, z)
        return z
    
    def self_consistent( self, T, Pk, ave_k ):
        '''Solves a self consistent equation using the G_1
        generating function as `u = G_1(u;T). Modifications 
        can be made where required.
        
        :param T: Transmissibility 
        :param Pk: degree distribution
        :param ave_k: average degree'''
        u = 0
        for i in range(1000):
            u = 1 - T + T*u
            u = self.G_1_generating_function(Pk, ave_k, self.return_arg, u)
        return u
        
    def do( self, params ):
        '''Runs the experiment. In this case we have an SIR model
        where a disease spreads over the network. 
        :param params: the experimental parameters'''
        kmean = params[self.AVERAGE_K]
        N = params[self.N]
        T = params[self.T]
        
        rc = dict()

        # create network
        g = self._network

        # find degree distribution
        Pk = self.degree_distribution(g)

        # find the average degree
        ave_k = self.average_degree(Pk)
        
        # compute final size of epidemic
        u = self.self_consistent(T, Pk, ave_k)     
        S_1 = self.super_critical_GC(T, Pk, ave_k, u)  
        
        rc['Pk'] = Pk 
        rc['ave_k'] = ave_k
        
        rc['u'] = u
        rc['S_1'] = S_1
        
        return rc