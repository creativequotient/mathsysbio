"""
This is the Main driver script that will handle the simulation of the system
"""
from BactSim.Bacteria import make_basic_bacteria
from BactSim.Simulator import Simulator, IntSimulator
from BactSim.FoodGenerators import FoodGenerator
from BactSim.FoodGenerators import StaticGenerator

import networkx as nx
import matplotlib.pyplot as plt

initial_bacteria = make_basic_bacteria(1)

food_source = StaticGenerator()
simulator = IntSimulator(food_generator=food_source,
                         initial_bacteria=[initial_bacteria])

print(initial_bacteria)

#print(initial_bacteria.get_edge("glucose", "transported_glucose"))

initial_bacteria = make_basic_bacteria(1)

#print(initial_bacteria.graph.nodes['atp']['amount'])
while True:
    simulator.progress()
    print(len(simulator.bacteria))
    print(simulator.bacteria[0])
    input("Continue?")
