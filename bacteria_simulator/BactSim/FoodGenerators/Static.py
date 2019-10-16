import math


default_food = {"glucose":1000,
                "sucrose":0,
                "lactose":0}

class StaticGenerator(object):

    def __init__(self, food = default_food):
        self.food = food
        self.generation = 0

    def getAvailable(self):
        return self.food  # returns a dictionary of amounts of each type of sugar
