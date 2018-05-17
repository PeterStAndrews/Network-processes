# stochastic heterogeneous mean field theory 
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

import epyc
import numpy as np

class STO( epyc.Experiment ):
    '''Stochastic simulation for a Markov process over an epyc experiment.
    Must specify the list of possible events and their effects on the system 
    via a subclass.'''
    
    MAX_TIME = 1000
    
    def initialisation( self, params ):
        '''Returns array of inititalised model. 
        :param params: the experimental parameters'''
        raise NotImplementedError 
       
    def after( self, states, rc):
        '''Function to be run after the simulation
        records the final populations at equilibrium or when 
        no events are left as a dict. 
        :param states: the final states
        :param rc: results dict'''
        raise NotImplementedError 
        
    def computeRates( self, params, states ):
        '''Computes the event rates and appends the event list.
        :param params: the experimental parameters
        :param states: the current state, array'''
        raise NotImplementedError 
        
    def transition_matrix( self ): 
        '''Returns transition matrix between states as array.'''
        raise NotImplementedError 
        
    def draw( self, t, params, states ):
        '''A single step of the algorithm, draws an event index and its time.
        Returns event rate index k and time t pair.
        :param t: current time
        :param params: experimental parameters
        :param states: current state'''
        e = None
        
        # compute event rates
        e_list = self.computeRates(params, states)
        
        # sum event rates
        sum_e = sum( e_list )
        
        # if event rate is None
        if sum_e == 0.0:
            return e, t
        
        else:
            # draw 1st random number
            u1 = np.random.uniform(0,1)
            
            # select event type e
            tot = u1 * sum_e
            e = 0
            k = e_list[e] 
            
            while k < tot:
                e += 1
                ev = e_list[e]
                k += ev 
    
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
        '''Performs the simulation until convergence through max time or 
        sum of event rates reaching zero (no more events left).
        :param parmas: experimental paramteres'''
        
        t = 0.
        rc = dict()
        states = self.initialisation(params)
        update_matrix = self.transition_matrix()
        
        while not self.at_equilibrium(t):
            
            # choose an event e and its time t
            e, t = self.draw(t, params, states)
            
            # check if any events left
            if e is None:
                break;
                
            # update model 
            states += update_matrix[e]
        
        self.after(states, rc)
        
        return rc