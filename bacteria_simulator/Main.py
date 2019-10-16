"""
This is the Main driver script that will handle the simulation of the system
"""
from BactSim.Bacteria import make_basic_bacteria
from BactSim.Simulator import Simulator
from BactSim.FoodGenerators import FoodGenerator

initial_bacteria = make_basic_bacteria(1)
food_source = FoodGenerator()
simulator = Simulator(food_generator=food_source,
                      initial_bacteria=initial_bacteria)
for i in range(5):
    simulator.progress()