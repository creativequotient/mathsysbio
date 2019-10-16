import networkx as nx
import copy
from BactSim.Evolution import Evolution

def make_basic_bacteria(id):
    def make_edge_cfg(weight, scale = 1, atp = 0, evo_sd = 0):
        return {
            'weight' : weight,
            'scale' : scale,
            'atp': atp,
            'make_evolution_cls' : lambda w : Evolution(w, evo_sd)
        }

    def make_all_edge_cfgs(transport, es, atp):
        return {
            'transport' : transport,
            'es_complex' : es,
            'atp' : atp
        }

    bac = Bacteria(id, min_atp = 30, initial_atp = 65)
    
    foods = ('glucose', 'sucrose', 'lactose')

    cfgs = {
        'glucose' : make_all_edge_cfgs(
            make_edge_cfg(0.8, evo_sd = 0.02, atp = 0.5),
            make_edge_cfg(0.8, evo_sd = 0.02, atp = 0.5),
            make_edge_cfg(0.8, evo_sd = 0.02, scale = 6)
        ),
        'sucrose': make_all_edge_cfgs(
            make_edge_cfg(0.000000000001, evo_sd = 0.001, atp = 0.5),
            make_edge_cfg(0.000000000001, evo_sd = 0.001, atp = 0.5),
            make_edge_cfg(0.000000000001, scale = 6, evo_sd = 0.001)
        ),
        'lactose' : make_all_edge_cfgs(
            make_edge_cfg(0.7, evo_sd = 0.02, atp = 0.5),
            make_edge_cfg(0.5, evo_sd = 0.02, atp = 0.5),
            make_edge_cfg(0.5, scale = 6, evo_sd = 0.02)
        )
    }
    # TODO: add interactions btw different transporters/enz
    # assumption : amount of food is limiting, enzymes/transporters are in excess
    for food in foods:
        # nodes
        transported = 'transported_' + food
        es_complex = 'enz_' + food + '_complex'

        bac.add_node(food)
        bac.add_node(transported)
        bac.add_node(es_complex)
        bac.add_edge(food, transported, **cfgs[food]['transport'])
        bac.add_edge(transported, es_complex, **cfgs[food]['es_complex'])
        bac.add_edge(es_complex, 'atp', **cfgs[food]['atp'])

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

    def add_edge(self, src, dest, weight, make_evolution_cls = lambda w: Evolution(w, 0), atp = 0, scale = 1, description = ''):
        if src not in self.graph or dest not in self.graph:
            raise ValueError('src or dest node has not been created')
        self.graph.add_edge(src, dest, weight=weight, atp_needed = atp, evolution = make_evolution_cls(weight), scale = scale, description = description)

    def is_alive(self):
        return self.graph.nodes['atp']['amount'] >= self.min_atp

    def survive(self, food, reset_food = True):
        """
        Determine if cell survives self generation given the quantities of food available

        Args:
            food : dict of food source to amount
        """
        for (food_src, amount) in food.items():
            self.graph.nodes[food_src]['amount'] = amount
        self.next_timestep()

        # reset food amount
        if reset_food:
            for food_src in food:
                self.graph.nodes[food_src]['amount'] = 0

        return self.is_alive()

    ## Increment timestep by 1, update graph
    def next_timestep(self):
        self.timestep += 1

        # update the graph
        increase_by = {}
        for (src, dest, data) in self.graph.edges.data():
            if dest not in increase_by:
                increase_by[dest] = 0
            dest_produced = data['weight'] * self.graph.nodes[src]['amount'] * data['scale']
            increase_by[dest] += dest_produced

            self.graph.nodes['atp']['amount'] -= data['atp_needed'] * dest_produced

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

    def edges_str(self):
        string = ''
        for src, dest, data in sorted(self.graph.edges.data()):
            string += f"{src} {dest}. weight: {data['weight']} scale: {data['scale']} atp needed: {data['atp_needed']} {'description:' + data['description'] if data['description'] else ''}\n"
        return string

    def nodes_str(self):
        string = ''
        for node_name in sorted(self.graph.nodes):
            string += f"{node_name}: {self.graph.nodes[node_name]}\n"
        return string

    def __str__(self):
        string = f'#### Bacteria {self.id}. Generation: {self.generation}. Age: {self.timestep}. ATP threshold: {self.min_atp}\n'
        string += '## NODES\n'
        string += self.nodes_str()
        string += '## EDGES\n'
        string += self.edges_str()
        # string += f"\nCurrent ATP: {self.graph.nodes['atp']['amount']}\n"
        # for node_name in sorted(self.graph.nodes):
            # if node_name == 'atp':
                # continue
            # string += f"{node_name}: {self.graph.nodes[node_name]['amount']}\n"
        return string

if __name__ == '__main__':
    bac1 = make_basic_bacteria(1)
    print(bac1)

    bac2 = bac1.divide(2)
    print(bac2)
    
    print(bac1 is bac2)
    print(bac1.graph is bac2.graph)
    print(bac1.graph.nodes['glucose'] is bac2.graph.nodes['glucose'])
    print()

    food = {'glucose' : 5, 'sucrose' : 10, 'lactose' : 3}
    print('Food', food)
    print(bac1)
    for i in range(4):
        bac1.survive(food)
        print(bac1)
