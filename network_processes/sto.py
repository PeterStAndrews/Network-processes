# stochastic heterogeneous mean field theory 
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
import epyc
import numpy as np
import networkx

class STO( NETWORK ):
    '''Stochastic simulation for a degree based mean field system over a network.
    Must specify the list of possible events and their effects on the system 
    via a subclass. The algorithm draws an event and a time for it to occur
    before it updates the states. It does this for each degree independently.
    To modify the experiment the user should overide the following:
    
    :func initialisation: returns an array of the initial states for the 
    kth system.
    
    :func transition_matrix: returns the transition matrix for each event type
    as an np.array.
    
    :func compute_rates: returns an array of event propensities for the current state
    of the process. Event rates can be zero, in fact the experiment has equilibrated 
    when all events have zero propensity.'''
    
    MAX_TIME = 5000

    def draw( self, t, params, state, k, ave_k, Pk ):
        '''A single step of the algorithm, draws an event index and its time.
        Returns event rate index e and time t pair.
        
        :param t: current time
        :param params: experimental parameters
        :param state: current state
        :param k: the degree
        :param ave_k: average degree
        :param Pk: degree distribution
        
        :returns int, float: event index, time'''
        e = None
        
        # compute event rates
        e_list = self.compute_rates(params, state, k, ave_k, Pk)
    
        # sum event rates
        sum_e = sum(e_list)
        
        # if event rate is zero or negative propensity break
        if sum_e == 0.0 or any(x < 0.0 for x in e_list):
            return e, t
        
        else:
            # draw 1st random number
            u1 = np.random.uniform(0,1)
            
            # select event type e
            tot = u1 * sum_e
            e = 0
            i = e_list[e] 
            
            while i < tot:
                e += 1
                ev = e_list[e]
                i += ev 
    
            # calculate the timestep delta
            u2 = np.random.uniform(0,1)
            dt = ( 1.0 / sum_e ) * np.log( 1.0 / u2 )
            
            # increment time
            t += dt
            return e, t
    
    def at_equilibrium( self, t ):
        '''Returns true if default end time reached.
        :param t: current time'''
        if t > self.MAX_TIME:
            return True

    def do( self, params ):
        '''Performs the simulation until convergence is reached either through max time out or 
        sum of event rates reaching zero (no more events left). We first compute network properties 
        and the transition matrix before iterating through the degrees and integrating each system.
        The results are then mapped to macro network values.
        
        :param parmas: experimental paramteres'''
        
        t = 0.
        rc = dict()
        rec = dict()
        
        # grab a copy of the network & compute Nk
        g = self._network
        
        # find the degree distribution
        Pk = self.degree_distribution(g)
        
        # find the average degree of the network
        ave_k = self.average_degree(Pk)
        
        # define the transition matrix
        update_matrix = self.transition_matrix()
        
        # initialise the macro system
        states = self.initialisation(params, g)
        
        rc['y0'] = states
        
        for k in states.keys():
            
            # initialise the kth system
            state = states[k]
            
            while not self.at_equilibrium(t):
                
                # choose an event e and its time t
                e, t = self.draw(t, params, state, k, ave_k, Pk)
                
                # check if any events left
                if e is None:
                    break;
                
                # update model 
                state += update_matrix[e]
                 
            # record final kth states
            rec[k] = state
        
        # sum over the k-systems to macro values
        rc['final_state'] = map(np.sum, zip(*rec.values()))
        
        return rc
    
    def initialisation( self, params, g ):
        '''Inititalisation of populations in model. The nodes 
        in the network are seeded with probability `pInfected`, 
        their degree is stored and the number of infecteds of 
        that degree is recorded. The initial state matrix is then 
        created and returned as a dictionary.
        
        :param params: the experimental parameters
        :param g: the network
        :returns dict of arrays: states'''
        pInfected = params['pInfected']
        Ik = {} # number of infected seed nodes of degree k
        Nk = {} # number of nodes of degree k
        
        for node in g.nodes_iter():
            k = g.degree(node)
            Nk[k] = Nk.get(k,0) + 1
            if pInfected > np.random.random_sample():
                Ik[k] = Ik.get(k,0) + 1
    
        # create initial state dict
        y0 = {}
        for k in Nk.keys():
            if k in Ik:
                y0[k] = np.array([Nk[k] - Ik[k], Ik[k], 0], dtype=int)
            else:
                # force it?
                Nk[k] -= 1
                Ik[k] = 1
                y0[k] = np.array([Nk[k], Ik[k], 0], dtype=int)
                
                # don't force it!
                # y0[k] = np.array([Nk[k], 0.0, 0.0]) 
                
        return y0
       
    def transition_matrix( self ): 
        '''Returns transition matrix between states. In this 
        case, we have: T = [[infection event],[recovery event]].
        :returns array: transition matrix'''
        return np.array([[-1, 1, 0],
                         [ 0, -1, 1] ], dtype=np.int)
        
    def compute_rates( self, params, state, k, ave_k, Pk ):
        '''Computes the event rates and appends the event list.
        
        :param params: the experimental parameters
        :param state: the current state, array
        :param k: the degree
        :param ave_k: average degree
        :param Pk: degree distribution
        
        :returns list: event list
        '''
        # unpack the parameters
        pInfect = params['pInfect']
        pRecover = params['pRecover']
        
        # unpack the states
        S, I, R = state
        
        # compute total number of nodes
        N = S + I + R
        
        def theta( ave_k, Pk, I_k ):
            '''Function makes theta(t) for a non-correlated network.
            
            :param ave_k: average degree
            :param Pk: degree distribution
            :param i: I_k value
            :returns float: theta(t)'''
            summation = 0
            for k in Pk.keys():
                summation += (k - 1) * Pk[k] * I_k
            return ( summation + 0.0 ) / ave_k 
        
        # create the events
        e1 = (k * pInfect * S * theta(ave_k, Pk, I) + 0.0) / N
        e2 = pRecover * I
        
        return [e1, e2] 
