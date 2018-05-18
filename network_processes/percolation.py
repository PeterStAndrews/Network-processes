import epyc 
from network_processes import *
import numpy as np
import networkx


class PERCOLATION( NETWORK ):
    '''Base for conducting experimental percolation results on
    networks. 
    
    :References:
    -------------
    .. [1] Newman. (2002) 'The spread of epidemic disease on networks'
    .. [2] Miller. (2007) 'Epidemic size and probability in populations 
       with heterogeneous infectivity and susceptibility'
    '''
    
    T = 'T'
    
    def do( self, params ):
        '''Here we perform a percolation experiment on the network g.
        We remove nodes with a probability (1 - T) and return the 
        largest connected component in the resulting network.
        :param params: experimental parameters'''
        
        rc = dict() 
        
        # grab a copy of the network
        g = self._network
        
        # unpack required parameters 
        T = params[self.T]
        N = params[self.N]
        
        # iterate over network nodes and remove with probability 1 - T_eff
        es = []
        for e in g.edges():
            if np.random.random() < (1 - T):
                # add edge e to ebunch
                es.append(e)
                    
        # remove edges in ebunch from network
        g.remove_edges_from(es)
        
        # find the giant connected component of the network
        gc = max(networkx.connected_component_subgraphs(g), key=len)
        
        # report the occupied fraction of the network
        rc['occupied_fraction'] = (len(gc) + 0.0) / N
        
        return rc