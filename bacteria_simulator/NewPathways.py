def func(t):  # Function takes in cell 'time' and acting edge
    x_0 = np.array([0.5, 0.1, 0.001, 0.5, 0.5, 0.5])  # Define initial conditions
    A = 1 / x_0 - 1
    return 1 / (A * np.exp(-t) + 1)  # Returns array of values for weight of edge


def invfunc(w):  # Function takes in edge weight and acting edge
    x_0 = np.array([0.5, 0.1, 0.001, 0.5, 0.5, 0.5])  # Define initial conditions
    return np.log(abs((1 / x_0 - 1) / (1 / w - 1)))  # Returns cell 'time'


def err(w, sd):  # Takes standard deviation of gaussian error as input
    return np.random.normal(scale=sd, size=len(w))


def NewPathways(w):  # Function takes 1x6 array [glu_dep,lac_dep,suc_dep,glu_eff,lac_eff,suc_eff]
    sd = 0.001  # Adjust the standard deviation of the Gaussian error
    print(func(invfunc(w) + err(w, sd)))
    return func(invfunc(w) + err(w, sd))
