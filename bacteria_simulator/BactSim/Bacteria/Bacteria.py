import networkx as nx
import copy
from BactSim.Evolution import Identity

def make_basic_bacteria(id):
    bac = Bacteria(id, min_atp = 5, initial_atp = 10)
    
    foods = ('glucose', 'sucrose', 'lactose')

    # TODO: tweak initial values, add interactions btw different transporters/enz
    # assumption : amount of food is limiting, enzymes/transporters are in excess
    for food in foods:
        bac.add_node(food)
        bac.add_node('transported_' + food)
        bac.add_node('enz_' + food + '_complex')
        bac.add_edge(food, 'transported_' + food, 0.5, atp = 1)
        bac.add_edge('transported_' + food, 'enz_' + food + '_complex', 0.5, atp = 1)
        bac.add_edge('enz_' + food + '_complex', 'atp', 0.5)

    return bac

class Bacteria(object):
    def __init__(self, id, min_atp, initial_atp):
        """
        Initalize single Bacterial cell with the given ATP threshold for survival (min_atp) 
        and initial amount of ATP (intitial_atp)
        Only the ATP node is created
        """
        self.graph = nx.DiGraph()
        self.graph.add_node('atp', amount=initial_atp)
        self.min_atp = min_atp
        self.initial_atp = initial_atp

        self.id = id
        self.timestep = 0
        self.generation = 0

    def add_node(self, name, initial_amount = 0, description = ''):
        self.graph.add_node(name, amount = initial_amount, descripton = '')

    def add_edge(self, src, dest, weight, make_evolution_cls = Identity, atp = 0, description = ''):
        if src not in self.graph or dest not in self.graph:
            raise ValueError('src or dest node has not been created')
        self.graph.add_edge(src, dest, weight=weight, atp_needed = atp, evolution = make_evolution_cls(weight), description = description)

    def survive(self, food):
        """
        Determine if cell survives self generation given the quantities of food available

        Args:
            food : dict of food source to amount
        """
        for (food_src, amount) in food.items():
            self.graph.nodes[food_src]['amount'] = amount
        self.next_timestep()

        # reset food amount
        # for food_src in food:
            # self.graph.nodes[food_src]['amount'] = 0

        return self.graph.nodes['atp']['amount'] >= self.min_atp

    ## Increment timestep by 1, update graph
    def next_timestep(self):
        self.timestep += 1

        # update the graph
        increase_by = {}
        for (tf, gene, data) in self.graph.edges.data():
            if gene not in increase_by:
                increase_by[gene] = 0
            gene_produced = data['weight'] * self.graph.nodes[tf]['amount']
            increase_by[gene] += gene_produced

            self.graph.nodes['atp']['amount'] -= data['atp_needed'] * gene_produced

        # increase the amounts for each node
        for (node, amount) in increase_by.items():
            if amount > 0:
                # it's not possible for a TF to actually decrease the amount of protein
                # only decrease its rate of production (?)
                self.graph.nodes[node]['amount'] += amount

    def clone(self, id):
        """Clone cell. Inherits age and amounts of nodes."""

        cloned_cell = copy.deepcopy(self)
        cloned_cell.id = id
        cloned_cell.generation += 1
        return cloned_cell

    def evolve(self):
        """Evolves (ie. possibly mutates) weights of all edges in this bacterium's graph"""

        for (src, dest, evolution) in self.graph.edges.data('evolution'):
            self.graph.edges[src, dest]['weight'] = evolution.getMutated()

    def divide(self, id):
        """Equivalent to this.clone(id) and using evolve() on the resulting daughter cell."""
        daughter = self.clone(id)
        daughter.evolve()
        return daughter

    def __str__(self):
        string = f'Bacteria {self.id}. Generation: {self.generation}. Age: {self.timestep}. ATP threshold: {self.min_atp}'
        string += f"\nCurrent ATP: {self.graph.nodes['atp']['amount']}\n"
        for node_name in sorted(self.graph.nodes):
            if node_name == 'atp':
                continue
            string += f"{node_name}: {self.graph.nodes[node_name]['amount']}\n"
        return string

if __name__ == '__main__':
    bac1 = make_basic_bacteria(1)
    print(bac1)
    print(bac1.graph.nodes.data())
    print(bac1.graph.edges.data())
    print()

    bac2 = bac1.divide(2)
    print(bac2)
    
    print(bac1 is bac2)
    print(bac1.graph is bac2.graph)
    print(bac1.graph.nodes['glucose'] is bac2.graph.nodes['glucose'])
    print()

    food = {'glucose' : 5, 'sucrose' : 10, 'lactose' : 3}
    print('Food', food)
    print(bac1.survive(food))
    print(bac1)
    for i in range(4):
        bac1.next_timestep()
        print(bac1)
