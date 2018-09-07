# Generating functions base class
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


class GFs( NETWORK ):
    '''Base for conducting epidemic analytical generating functions on
    networks. This class defines common methods that will be used frequently. 
    Sub-classes will define the details of the experimental sections and will 
    overide `do` in all but the basic scenarios. Currently, the generating 
    functions are designed for the configuration model; however, it is clear 
    how to customise these when required.
    
    :References:
    -------------
    .. [1] Newman. (2002) 'The spread of epidemic disease on networks'
    .. [2] Miller. (2007) 'Epidemic size and probability in populations 
       with heterogeneous infectivity and susceptibility'
    .. [3] Funk & Jansen (2010) `Interacting epidemics on overlay networks`
    '''     
    T = 'T' # transmissibility 
    
    def __init__(self):
        super(GFs, self).__init__()
        
    def G_0_generating_function( self, Pk, ave_k, x ):
        r'''Computes a `G_0` generating function as 
        
         .. math::
         
             \sum_{k=0}^\infty p_k(x)^k
        
        :param Pk: degree distribution
        :param ave_k: average degree
        :param x: argument function'''
        G_0 = 0
        for k in Pk.keys():
            G_0 += Pk[k]*x**k
        return G_0
    
    def G_1_generating_function( self, Pk, ave_k, x ):
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
                summation += k*Pk[k]*x**(k-1)
        return (summation + 0.0)/ave_k
    
    def do( self, params ):
        '''Runs the experiment.  '''
        T1 = params[self.T]
        
        rc = dict()
        g = self._network
        Pk = self.degree_distribution(g)
        ave_k = self.average_degree(Pk)
        
        # first disease 
        u = 0.5
        for i in range(0,2000):
            u = 1 - self.G_1_generating_function(Pk, ave_k, 1-u*T1)
        S1 = 1 - self.G_0_generating_function(Pk, ave_k, 1-u*T1)
    
        rc['S1'] = S1
        rc['Pk'] = Pk 
        rc['ave_k'] = ave_k
        
        return rc
