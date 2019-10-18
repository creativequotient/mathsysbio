"""
This is the Main driver script that will handle the simulation of the system
"""
from BactSim.Bacteria import make_basic_bacteria
from BactSim.Simulator import Simulator, IntSimulator
from BactSim.FoodGenerators import FoodGenerator
from BactSim.FoodGenerators import StaticGenerator

initial_bacteria = make_basic_bacteria(1)
food_source = StaticGenerator()
simulator = IntSimulator(food_generator=food_source,
                         initial_bacteria=[initial_bacteria])
for i in range(20):
    simulator.progress()
    print(len(simulator.bacteria))
