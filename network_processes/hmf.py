# Heterogeneous mean field theory base class
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
from scipy.integrate import ode



class HMF( epyc.Experiment ):
    '''This class will generate a network and from it obtain the 
    degree distribution. This will be used to integrate a mean 
    field SIR system according to the heterogeneous mean field 
    theory.
    
    :References:
    -------------
    .. [1] R. Pastor-Satorras and A. Vespignani. 'Epidemic spreading 
           in scale-free networks', Phys. Rev. Lett., vol. 86, pp. 3200-3203, 2001. '''
    
    N = 'N' # order of the network
    AVERAGE_K = 'kmean' # average degree s
    
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

    def do( self, params ):
        '''runs the experiment.'''
        rc = dict()
        
        pInfect = params['pInfect']
        pRecover = params['pRecover']
        pInfected = params['pInfected']
        
        S = dict()
        I = dict()
        R = dict()
        
        # grab a copy of the network
        g = self._network
        
        # find the degree distribution
        Pk = self.degree_distribution(g)
        
        # find the average degree of the network
        ave_k = self.average_degree(Pk)
        
        def theta( ave_k, Pk, I_k ):
            '''Function makes theta(t) for a non-correlated network.
            
            :param ave_k: average degree
            :param Pk: degree distribution
            :param i: I_k value
            :returns: theta(t)'''
            summation = 0
            for k in Pk.keys():
                summation += (k - 1) * Pk[k] * I_k
            return ( summation + 0.0 ) / ave_k 
    
        def make_sir( t, y, pInfect, pRecover, k, ave_k, Pk ):
            '''Return functions for changes in the susceptible, infectious, and recovered
            sub-populations for particular rates of recovery and infection.
            
            :param pInfect: rate of infection
            :param pRecover: rate of recovery
            :returns: change functions'''
            
            S, I, R = y
            
            dS = - k * pInfect * S * theta(ave_k, Pk, I)
            dI = k * pInfect * S * theta(ave_k, Pk, I) - pRecover * I
            dR = pRecover * I
            
            return [dS, dI, dR]
        
        for k in Pk.keys():
            
            # initial conditions vector for kth system
            y0 = [1 - pInfected, pInfected, 0]
            t0 = 0
            
            t1 = 150 
            dt = 1
            
            r = ode(make_sir).set_integrator('dopri5', method = 'adams')
            r.set_initial_value(y0, t0).set_f_params(pInfect, pRecover, k, ave_k, Pk)       
            
            while r.successful() and r.t < t1:
                r.integrate(r.t+dt)
            
            # report final results for kth system
            S[k] = r.y[0] * Pk[k]
            I[k] = r.y[1] * Pk[k]
            R[k] = r.y[2] * Pk[k]
        
        rc['S_final'] = sum(S.values())
        rc['I_final'] = sum(I.values())
        rc['R_final'] = sum(R.values())
        
        return rc