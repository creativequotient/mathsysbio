class Bacteria(object):

    ## Initalize single Bacterial cell with the the given weights and input
    def __init__(self, ID, generation, edges, mutations):
        this.edges = edges
        this.ID = ID
        this.generation = generation
        this.mutations = mutations
        '''
        Maybe come up with some way to couple the objects in mutations (which correspond to the mutation functions for each edge) ;)
        '''

    ## Determine if cell survives this generation given the quanitities of food available
    def survive(self, food):
        return True

    ## Clone cell
    def clone(self):
        return 0

    ## Evolve daughter cell
    def evolve(self):
        return 0
