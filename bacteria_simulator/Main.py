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


def glucose_stats(bacteria_pop):
    output = [0,0,0]
    for bac in bacteria_pop:
        output[0] += bac.get_edge("glucose", "transported_glucose")["weight"]
        output[1] += bac.get_edge("transported_glucose", "enz_glucose_complex")["weight"]
        output[2] += bac.get_edge("enz_glucose_complex", "atp")["weight"]
    output = list(map(lambda x: x / len(bacteria_pop), output))
    print("Glucose -> Transported Glucose {0:.5f} Transported_glucose -> Enz_Glucose_Complex {1:.5f} Enz_Glucose_Complex -> ATP {2:.5f}"
          .format(output[0][0], output[1][0], output[2][0]))

def lactose_stats(bacteria_pop):
    output = [0,0,0]
    for bac in bacteria_pop:
        output[0] += bac.get_edge("lactose", "transported_lactose")["weight"]
        output[1] += bac.get_edge("transported_lactose", "enz_lactose_complex")["weight"]
        output[2] += bac.get_edge("enz_lactose_complex", "atp")["weight"]
    output = list(map(lambda x: x / len(bacteria_pop), output))
    print("Lactose -> Transported Lactose {0:.5f} Transported_lactose -> Enz_Lactose_Complex {1:.5f} Enz_Lactose_Complex -> ATP {2:.5f}"
          .format(output[0][0], output[1][0], output[2][0]))

def sucrose_stats(bacteria_pop):
    output = [0,0,0]
    for bac in bacteria_pop:
        output[0] += bac.get_edge("sucrose", "transported_sucrose")["weight"]
        output[1] += bac.get_edge("transported_sucrose", "enz_sucrose_complex")["weight"]
        output[2] += bac.get_edge("enz_sucrose_complex", "atp")["weight"]
    output = list(map(lambda x: x / len(bacteria_pop), output))
    print("Sucrose -> Transported Sucrose {0:.5f} Transported_sucrose -> Enz_Sucrose_Complex {1:.5f} Enz_Sucrose_Complex -> ATP {2:.5f}"
          .format(output[0][0], output[1][0], output[2][0]))

#print(initial_bacteria)

for i in range(10000):
    simulator.progress()
    if i % 10 == 0:
        print("Food available: {}".format(food_source.food))
        print("Generation: {}".format(i))
        print("Population size: {}".format(len(simulator.bacteria)))
        glucose_stats(simulator.bacteria)
        lactose_stats(simulator.bacteria)
        sucrose_stats(simulator.bacteria)
        input("Continue")
