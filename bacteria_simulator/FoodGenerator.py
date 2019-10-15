import math

class FoodGenerator(object):

    def __init__(self):
        self.food = {"glucose": 0,
                     "lactose": 0,
                     "sucrose": 0}
        self.generation = 0

    def getAvailable(self):
        """
        Get available food at current generation
        :return: Dictionary of food, value pairs
        """
        self.food["glucose"] = 3 * math.sin((self.generation - 3)/100) + 3  # to generate the amount of glucose in generation self.generation
        self.food["lactose"] = 3 * math.sin((self.generation - 5)/100) + 3  # to generate the amount of lactose in generation self.generation
        self.food["sucrose"] = 3 * math.sin((self.generation - 7)/100) + 3  # to generate the amount of sucrose in generation self.generation
        self.generation += 1
        return self.food  # returns a dictionary of amounts of each type of sugar


if __name__ == "__main__":
    fg = FoodGenerator()
    for i in range(200):
        print(fg.getAvailable())