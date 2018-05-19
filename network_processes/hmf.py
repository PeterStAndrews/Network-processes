# Degree based heterogeneous mean field theory RK4 integration
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

from network_processes import *
import networkx
import epyc
from scipy.integrate import ode
import numpy as np

class HMF( NETWORK ):
    '''This class will integrate a degree based mean field system using an RK4
    integration scheme. The network is created in the base class `NETWORK` The user 
    should modify `model`, `initialisation` and `param_vector` to tailor this to their needs. 
    
    :func model: returns the updated differential equations dx/dt as a np.array. 
    :func initialisation: returns the initial state {X_0} for the kth model as an
    an np.array. 
    :func param_vector: returns the parameters required to update the model using `model`.
    
    :References:
    -------------
    .. [1] R. Pastor-Satorras and A. Vespignani. 'Epidemic spreading 
           in scale-free networks', Phys. Rev. Lett., vol. 86, pp. 3200-3203, 2001. '''
    
    def model( self, t, y, pInfect, pRecover, k, ave_k, Pk ):
        '''Return functions for changes in states for the kth model.
        :param pInfect: rate of infection
        :param pRecover: rate of recovery
        :returns: change functions'''
            
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
            
        # unpack current state
        S, I, R = y
            
        # define the update equations 
        dS = - k * pInfect * S * theta(ave_k, Pk, I)
        dI = k * pInfect * S * theta(ave_k, Pk, I) - pRecover * I
        dR = pRecover * I
            
        return np.array([dS, dI, dR])
    
    def initialisation( self, params):
        '''Initialises the kth state.
        :param params: the experimental parameters
        :returns array: the initial kth state'''
        pInfected = params['pInfected']
        return np.array([1 - pInfected, pInfected, 0], dtype=float)
        
    def param_vector( self, params, k, ave_k, Pk):
        '''A vector of parameters required for the integrator.
        :param params: the experimental parameters
        :param k: the degree
        :param ave_k: the average degree
        :param Pk: the degree distribution'''
        pInfect = params['pInfect']
        pRecover = params['pRecover']
        return pInfect, pRecover, k, ave_k, Pk
    
    def do( self, params ):
        '''runs the experiment.'''
        rc = dict()
        states = dict()
        
        # grab a copy of the network
        g = self._network
        
        # find the degree distribution
        Pk = self.degree_distribution(g)
        
        # find the average degree of the network
        ave_k = self.average_degree(Pk)

        for k in Pk.keys():
            
            # initial conditions vector for kth system
            y0 = self.initialisation(params)
            
            vec = self.param_vector(params, k, ave_k, Pk)
            t0 = 0
            
            t1 = 150 
            dt = 1
            
            r = ode(self.model).set_integrator('dopri5', method = 'adams')
            r.set_initial_value(y0, t0).set_f_params(*vec)       
            
            while r.successful() and r.t < t1:
                r.integrate(r.t+dt)
            
            # report final results for kth system
            states[k] = r.y * Pk[k]
        
        rc['final_state'] = sum(states.values())
        
        return rc
