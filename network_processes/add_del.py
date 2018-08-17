
# Addition-deletion process for an evolving network.
#
# Copyright (C) 2018 Peter Mann
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

from network_processes import *
import epyc
from scipy.integrate import odeint
import numpy as np
from collections import defaultdict

class addition_deletion(NETWORK):
    '''Integrates a dim(k_max,n) system of ODEs using
    `scipy.integrate.odeint`. To customise this class, 
    modify the `initialisation` method, update the rate 
    equation calculation and how the results are reported.
    '''
    def __init__(self):
        super(addition_deletion, self).__init__()
        
    def initialisation( self, G, k_max, n):
        '''Creates a dict (dim(k)) of lists (dim(n)) of total dim(k_max,n) 
        and populates with zeros. The first column is populated with the degree
        distribution for the network G and the dictionary values are flattened.
        :param G: the network
        :param k_max: max degree cut-off
        :param n: class dimension
        :returns p0: flattened vector of dim(k_max,n)'''
        # create dict of zeros
        p0 = defaultdict(list)
        for k in range(0, k_max):
            p0[k] = [0]*n
        
        # populate first column with p[k] values
        N = G.order()
        inv_N = 1./N
        for node in G.nodes_iter():
            k = G.degree(node)
            p0[k][0] += inv_N
        
        # flatten the dict values into an array
        return np.array([item for sublist in p0.values() for item in sublist]) 
    
    def dpdt( self, p, t, k_max, n, c ):
        '''Computes the rate equation for an addition-
        deletion network [1].
        
        :param y: flattened state vector
        :param t: current time
        :param k_max: the max cut-off degree
        :param n: class dimension
        :param c: average degree of new node
        
        :returns dy: flattened list of change functions
        
        :References:
        -------------
        .. [1] C. Moore, G. Ghoshal, and M. E. J. Newman, 
           “Exact solutions for models of evolving networks 
           with addition and deletion of nodes,” Phys. Rev. E, 
           vol. 74, p. 036121, Sep 2006.
        '''
        # then reshape it to the original dimension (k_max*n)
        p = p.reshape(k_max, n) 
    
        # a fresh dict to store dp_k values
        dp = defaultdict(list)
        for k in range(0, k_max):
            dp[k] = [0]*n   
        
        def _phi_k( c, k, POISSON, DELTA):
            '''Evaluates phi_k for arg set.
            :param c: the average degree of an oncoming node
            :param k: the degree
            :param POISSON: Bool flag for Poisson distribution
            :param DELTA: Bool flag for delta function, (c,k)'''
            if POISSON:
                return (np.exp(-c)*c**k)/np.math.factorial(k)
            
            if DELTA:
                if k is c:
                    return 1.0
                else:
                    return 0.0
            
        # iterate over ks
        for k, value in enumerate(p):
            
            # determine phi_k
            phi_k = _phi_k(c, k, self._POISSON, self._DELTA)
            
            # rate equations for p_k
            if k is 0:
                # equation for k=0 case
                dp[k] = -p[k] - k*p[k] + (k+1)*p[k+1] - c*p[k] + phi_k
                
            if k > 0 and k + 1 < len(p):
                # equation general case
                dp[k] = -p[k] - k*p[k] + (k+1)*p[k+1] + c*p[k-1] - c*p[k] + phi_k
            
            if k is len(p)-1:
                # equation for k=k_max case
                dp[k] = - p[k] - k*p[k] + c*p[k-1] - c*p[k] + phi_k
            
        return [item for sublist in dp.values() for item in sublist]
    
    def do( self, params ):
        '''Performs the experiment. Un-packs the 
        parameters and creates the network. The initialisation
        vector is populated and the system is integrated.
        :param params: experimental parameters'''
        # dict to store results 
        rc = dict()
        
        # un-pack the parameters
        t = params['time']
        n = params['class_dimension']
        N = params['N']
        k_max = params['k_max']
        c = params['kmean']
        self._POISSON = params['poisson']
        self._DELTA = params['delta']
        
        # grab a copy of the network
        G = self._network
        
        # initialise the sytem
        p0 = self.initialisation(G, k_max, n)
        
        # integrate the system
        soln = odeint(self.dpdt, p0, t, args=(k_max, n, c))
        
        # report the results
        rc['sol'] = soln[-1]
        
        return rc
    
