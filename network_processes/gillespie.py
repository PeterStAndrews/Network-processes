import epyc
import numpy as np
import networkx

class Gillespie( epyc.Experiment ):
    '''Stochastic simulation for a Markov process over an epyc experiment.
    Must specify the list of possible events and their effects on the system 
    via a subclass.'''
    
    MAX_TIME = 5000
    
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
        :param g: the network
        :returns dict: degree distribution'''
        Pk = {}
        order = g.order()
        inv_order = 1./order
        for node in g.nodes_iter():
            k = g.degree(node)
            Pk[k] = Pk.get(k,0) + inv_order
        return Pk

    def average_degree( self, Pk ):
        '''Returns the average degree of the degree distribution Pk.
        :param Pk: degree distribution
        :returns float: average degree'''
        ave_k = 0
        for k in Pk.keys():
            ave_k += k*Pk[k]
        return ave_k
    
    def degree_composition( self, g ):
        '''Populates N[k] dict containing number of nodes 
        of degree k.
        :param g: the network
        :returns dict: degree composition'''
        Nk = {}
        for node in g.nodes_iter():
            k = g.degree(node)
            Nk[k] = Nk.get(k,0) + 1
        return Nk
    
    def initialisation( self, params, N ):
        '''Returns array of inititalised model. 
        :param params: the experimental parameters
        :params N: number of nodes
        :returns array: initial states'''
        raise NotImplementedError 
       
    def after( self, states, rc):
        '''Function to be run after the simulation
        records the final populations at equilibrium or when 
        no events are left as a dict. 
        
        :param states: the final states
        :param rc: results dict'''
        raise NotImplementedError 
        
    def computeRates( self, params, states, k, ave_k, Pk ):
        '''Computes the event propensities and appends the event list.
        
        :param params: the experimental parameters
        :param states: the current state, array
        :param k: the degree
        :param ave_k: average degree
        :param Pk: degree distribution
        :returns list: event list'''
        raise NotImplementedError 
        
    def transition_matrix( self ): 
        '''Returns transition matrix between states as array.
        :returns array: transition matrix'''
        raise NotImplementedError 
        
    def draw( self, t, params, states, k, ave_k, Pk ):
        '''A single step of the algorithm, draws an event index and its time.
        Returns event rate index e and time t pair.
        
        :param t: current time
        :param params: experimental parameters
        :param states: current state
        :param k: the degree
        :param ave_k: average degree
        :param Pk: degree distribution
        
        :returns int, float: event index, time'''
        e = None
        
        # compute event rates
        e_list = self.computeRates(params, states, k, ave_k, Pk)
    
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
        Nk = self.degree_composition(g)
        
        # find the degree distribution
        Pk = self.degree_distribution(g)
        
        # find the average degree of the network
        ave_k = self.average_degree(Pk)
        
        # define the transition matrix
        update_matrix = self.transition_matrix()
        
        for k in Nk.keys():
            
            # initialise the kth system
            states = self.initialisation(params, Nk[k])
            
            while not self.at_equilibrium(t):
                
                # choose an event e and its time t
                e, t = self.draw(t, params, states, k, ave_k, Pk)
                
                # check if any events left
                if e is None:
                    break;
                
                # update model 
                states += update_matrix[e]
                 
            # record final kth states
            rec[k] = states
        
        # sum over the k-systems to macro values
        rc['final_states'] = map(np.sum, zip(*rec.values()))
        rc['Nk'] = Nk
        
        return rc