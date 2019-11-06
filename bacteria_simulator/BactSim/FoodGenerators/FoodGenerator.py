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
        rad = self.generation * math.pi / 180 / 0.2
        self.food["glucose"] = 500 * math.cos(rad) + 500 #+ 500  # to generate the amount of glucose in generation self.generation
        self.food["lactose"] = 500 * math.cos(0.9 * rad) + 500 #+ 500 to generate the amount of lactose in generation self.generation
        self.food["sucrose"] = 500 * math.sin(rad + 1.5 * math.pi) + 500 # to generate the amount of sucrose in generation self.generation
        self.generation += 1
        return self.food  # returns a dictionary of amounts of each type of sugar


if __name__ == "__main__":
    fg = FoodGenerator()
    for i in range(10000):
        if i % 1 == 0:
            print(fg.getAvailable())
