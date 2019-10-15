import numpy as np

class Evolution(object):

    def __init__(self, initial, sd):
        """
        Initialize the Evolution class, an object to handle mutation of weights
        :param initial: initial weight of edge at generation/time 0
        :param sd: standard deviation of errors
        """
        self.initial = initial
        self.sd = sd
        self.weight = initial

    def evolve(self, time):
        """
        Calculate weight given time
        :param time: time
        :return: weight
        """
        A = 1 / self.initial - 1
        return 1 / (A * np.exp(-time) + 1)

    def inverse(self, weight):
        """
        Returns corresponding time given weight
        :param weight: weight
        :return: corresponding time to weight
        """
        return np.log(abs((1 / self.initial - 1) / (1 / weight - 1)))

    def error(self):
        """
        Generate error given sd
        :return: error within self.sd
        """
        return np.random.normal(scale=self.sd, size=1)

    def getMutated(self):
        """
        Generate next mutated weight
        :return: new weight with slight mutation from the previous weight
        """
        self.weight = self.evolve(self.inverse(self.weight) + self.error())
        if (self.weight > 1):
            self.weight = 1
        return self.weight

if __name__ == "__main__":
    evo = Evolution(0.5, 0.001) # Initialize Evolution object with initial weight and SD
    for i in range(50):
        print(evo.getMutated()) # Call .getMutated() to get next weight