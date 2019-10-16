class Identity(object):

    def __init__(self, initial):
        """
        Evolution class for edges we don't want to allow changes in weight. Useful for testing
        :param initial: Initial weight of edge
        """
        self.weight = initial

    def getMutated(self):
        return self.weight