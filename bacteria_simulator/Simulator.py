class Simulator(object):

    def __init__(self, food_generator, initial_bacteria):
        """
        Initialize Simulator class
        :param food_generator: FoodGenerator object to output food available at each generation
        :param initial_bacteria: Initial bacteria cell at generation 0
        """
        self.food_generator = food_generator
        self.bacteria = [initial_bacteria]
        self.total_population = 1

    def progress(self):
        """
        Move forward by 1 generation/time-step
        :return: Update self.bacteria with new population at next time-step
        """
        new_population = self.replicate()
        available_food = self.calculate_food(len(new_population))
        self.bacteria = list(filter(lambda x: x.survive(available_food)))

    def replicate(self):
        """
        All members in current population replicate
        :return: List containing all new bacteria cells (cells from previous generation + progeny)
        """
        new_population = []
        for bacteria in self.bacteria:
            new_population.append(bacteria)
            self.total_population += 1
            daughter_cell = bacteria.clone().evolve()
            new_population.append(daughter_cell)
        return new_population

    def calculate_food(self, population_size):
        """
        Determine food available to each cell in population
        :param population_size: Size of current population
        :return: Dictionary containing food available to each cell
        """
        available_food = {}
        for food, quantity in list(self.food_generator.getAvailable().items()):
            available_food[food] = quantity / population_size
        return available_food

