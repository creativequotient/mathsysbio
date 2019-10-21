import random

class IntSimulator(object):
    """IntSimulator only allocates integer food values to bacteria as evenly as possible.
    All bacteria are first allocated amount // number of bacteria units of food.
    The remainder is then randomly allocated (at most 1 extra unit of food per bacteria).
    """

    def __init__(self, food_generator, initial_bacteria):
        """
        Initialize Simulator class
        :param food_generator: FoodGenerator object to output food available at each generation
        :param initial_bacteria: list of initial bacteria
        """
        self.food_generator = food_generator
        self.bacteria = initial_bacteria
        self.total_population = len(self.bacteria)

    def progress(self):
        """
        Move forward by 1 generation/time-step
        :returns: Update self.bacteria with new population at next time-step
        """
        new_population = self.replicate()
        food_alloc = self.food_allocation(len(new_population))
        self.bacteria = []
        for i, bac in enumerate(new_population):
            current_food = {}
            for food, lst in food_alloc.items():
                current_food[food] = lst[i]
            if bac.survive(current_food):
                self.bacteria.append(bac)

    def replicate(self):
        """
        All members in current population replicate
        :returns: List containing all new bacteria cells (cells from previous generation + progeny)
        """
        new_population = []
        for bacteria in self.bacteria:
            if not bacteria.can_reproduce():
                continue
            new_population.append(bacteria)
            self.total_population += 1
            daughter_cell = bacteria.clone(self.total_population)
            daughter_cell.evolve()
            new_population.append(daughter_cell)
        return new_population

    def food_allocation(self, population_size):
        """Generates food allocation scheme
        "param population_size: number of bacteria to allocate food to
        :returns: dict of sugar to list of food allocated to each bacterium"""
        available_food = {}
        for food, quantity in self.food_generator.getAvailable().items():
            min_qty = quantity // population_size
            lst = [min_qty] * population_size
            get_extra = int(quantity % population_size)
            for i in range(get_extra):
                lst[i] += 1
            random.shuffle(lst)
            available_food[food] = lst
        return available_food

# to run this file as a script:
# run BactSim.Simuator.IntSimulator from the bacteria_simulator directory
# ie. from this directory: cd ... && python3 -m BactSim.Simuator.IntSimulator
if __name__ == '__main__':
    from BactSim.FoodGenerators import StaticGenerator
    gen = StaticGenerator()
    simulator = IntSimulator(gen, [])
    print(gen.getAvailable().items())
    print(simulator.food_allocation(1))
    print(simulator.food_allocation(2))
    print(simulator.food_allocation(3))
