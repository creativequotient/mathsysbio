import random
from multiprocessing import Pool
import functools

class IntSimulator(object):
    """IntSimulator only allocates integer food values to bacteria as evenly as possible.
    All bacteria are first allocated amount // number of bacteria units of food.
    The remainder is then randomly allocated (at most 1 extra unit of food per bacteria).
    """

    def __init__(self, food_generator, initial_bacteria, food_unit = 10):
        """
        Initialize Simulator class
        :param food_generator: FoodGenerator object to output food available at each generation
        :param initial_bacteria: list of initial bacteria
        :param food_unit: allocate food in mutiples of this number (default: multiples of 10)
        """
        self.food_generator = food_generator
        self.bacteria = initial_bacteria
        self.total_population = len(self.bacteria)
        self.food_unit = food_unit

    def progress(self):
        """
        Move forward by 1 generation/time-step
        :returns: Update self.bacteria with new population at next time-step
        """
        if (len(self.bacteria) < 2000):
            new_population = self.replicate()
        else:
            new_population = replicate_multicore(self.bacteria)
        food_alloc = self.food_allocation(len(new_population))
        self.bacteria = []
        for i, bac in enumerate(new_population):
            current_food = {}
            for food, lst in food_alloc.items():
                current_food[food] = lst[i]
            #print(current_food)
            # print(bac)
            # print(bac.survive(current_food))
            if bac.survive(current_food):
                self.bacteria.append(bac)

    def replicate(self):
        """
        All members in current population replicate
        :returns: List containing all new bacteria cells (cells from previous generation + progeny)
        """
        new_population = []
        for bacteria in self.bacteria:
            #print(bacteria)
            if not bacteria.can_reproduce():
                continue
            new_population.append(bacteria)
            self.total_population += 1
            new_population.append(bacteria.divide(self.total_population))
        return new_population


    def food_allocation(self, population_size):
        """Generates food allocation scheme
        "param population_size: number of bacteria to allocate food to
        :returns: dict of sugar to list of food allocated to each bacterium"""
        food_unit = self.food_unit
        available_food = {}
        for food, quantity in self.food_generator.getAvailable().items():
            quantity /= food_unit
            min_qty = quantity // population_size
            lst = [food_unit * min_qty] * population_size
            get_extra = int(quantity % population_size)
            for i in range(get_extra):
                lst[i] += food_unit
            random.shuffle(lst)
            available_food[food] = lst
        #print(available_food)
        return available_food

def func(bacteria):
    if not bacteria.can_reproduce():
        return [bacteria]
    else:
        return [bacteria, bacteria.divide(0)]

def replicate_multicore(bacteria_pop, cores = 8):
    pool = Pool(processes=cores)
    output = pool.map(func, bacteria_pop)
    output = functools.reduce(lambda a, b: a + b, output)
    return output

# to run this file as a script:
# run BactSim.Simuator.IntSimulator from the bacteria_simulator directory
# ie. from this directory: cd ... && python3 -m BactSim.Simuator.IntSimulator
if __name__ == '__main__':
    from BactSim.FoodGenerators import StaticGenerator
    gen = StaticGenerator(food = {'glucose':139,'lactose':100})
    simulator = IntSimulator(gen, [], 10)
    print(gen.getAvailable().items())
    print(simulator.food_allocation(1))
    print(simulator.food_allocation(2))
    print(simulator.food_allocation(3))
    print(simulator.food_allocation(4))
