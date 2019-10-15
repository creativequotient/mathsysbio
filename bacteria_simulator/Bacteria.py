import networkx as nx
import copy


class Bacteria(object):
    def __init__(this, ID, generation, min_atp, initial_atp):
        """
        Initalize single Bacterial cell with the given ATP threshold for survival (min_atp) 
        and initial amount of ATP (intitial_atp)
        Only the ATP node is created

        Use add_product_node(), add_node() and add_edge() to add other nodes and edges.
        """
        this.graph = nx.DiGraph()
        this.graph.add_node('atp', amount=initial_atp)
        this.min_atp = min_atp
        this.initial_atp = initial_atp

        this.ID = ID
        this.timestep = 0
        this.generation = generation

    ## Functions to add nodes/edges to the graph
    def add_product_node(this, name, initial_amount, atp_needed):
        """
        Add a node a gene product and hence requires ATP to be synthesized

        Args:
            initial_amount (float) : initial amount of this gene product 
            atp_needed (float) : amount of ATP needed to produce each unit of this gene product
                (would it be better to associate the atp_needed with the pathway (ie. edge) 
                instead??? as different pathways may have different ATP requirements)
        """
        this.graph.add_node(name, amount=initial_amount, atp_fn=lambda produced: produced * atp_needed)

    ## Add non gene product nodes that don't need ATP (ie. food source nodes)
    def add_node(this, name, initial_amount):
        this.graph.add_node(name, amount=initial_amount)
        # print(this.graph.nodes[name]['amount'])

    def add_edge(this, tf, gene, weight, evolve_fn):
        if tf not in this.graph or gene not in this.graph:
            raise ValueError('tf or gene node has not been created')
        this.graph.add_edge(tf, gene, weight=weight, evolve_fn=evolve_fn)
        # print(this.graph.edges.data())

    ## Determine if cell survives this generation given the quantities of food available
    ## food is a dict of food_source : amount
    def survive(this, food):
        for (food_src, amount) in food.items():
            this.graph.nodes[food_src] = amount
            # this.graph.add_node(food_src,amount=amount)
        this._next_timestep()  # update the graph

        for food_src in food:
            this.graph.nodes[food_src]['amount'] = 0

        return this.graph.nodes['atp']['amount'] >= this.min_atp

    ## Increment timestep by 1, update graph
    def _next_timestep(this):
        this.timestep += 1

        # update the graph
        increase_by = {}
        for (tf, gene, weight) in this.graph.edges.data('weight'):
            if gene not in increase_by:
                increase_by[gene] = 0
            gene_produced = weight * this.graph.nodes[tf]['amount']
            # print("gene produced =" + str(gene_produced))
            increase_by[gene] += gene_produced

            # decrease ATP if gene product was produced
            if gene_produced > 0 and this._is_product(gene):
                this.graph.nodes['atp']['amount'] -= this.graph.nodes[gene]['atp_fn'](gene_produced)

        # increase the amounts for each node
        for (node, amount) in increase_by.items():
            if amount > 0:
                # it's not possible for a TF to actually decrease the amount of protein
                # only decrease its rate of production (?)
                this.graph.nodes[node]['amount'] += amount
        # print(increase_by)

    # it's a gene product (ptn/rna) if it has an atp function
    def _is_product(this, gene):
        return 'atp_fn' in this.graph.nodes[gene]

    ## Clone cell (inherits age and amounts of proteins)
    def clone(this):
        cloned_cell = copy.copy(this)
        cloned_cell.ID += 1
        return cloned_cell

    ## Evolve all edge in this cell
    def evolve(this):
        for (src, dest, evolve_fn) in this.graph.edges.data('evolve_fn'):
            this.graph.edges[src, dest]['weight'] = evolve_fn(this.graph.edges[src, dest]['weight'])

    def __str__(this):
        string = f'Bacteria {this.ID}. Generation: {this.generation}. Age: {this.timestep}. ATP threshold: {this.min_atp}'
        string += f"\nCurrent ATP: {this.graph.nodes['atp']['amount']}\n"
        for node_name in sorted(this.graph.nodes):
            if node_name == 'atp':
                continue
            string += f"{node_name}: {this.graph.nodes[node_name]['amount']}\n"
        return string
