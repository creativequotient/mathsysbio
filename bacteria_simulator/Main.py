"""
This is the Main driver script that will handle the simulation of the system
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
import time

from BactSim.Bacteria import make_basic_bacteria
from BactSim.Simulator import Simulator, IntSimulator
from BactSim.FoodGenerators import FoodGenerator
from BactSim.FoodGenerators import StaticGenerator

initial_bacteria = make_basic_bacteria(1)
food_source = StaticGenerator()
simulator = IntSimulator(food_generator=food_source,
                         initial_bacteria=[initial_bacteria],
                         food_unit=1)

def sugar_stats(bacteria_pop, sugar):
    """
    Get edge weights for each sugar pathway in the population
    :param bacteria_pop: Bacteria population
    :param sugar: Sugar/Food source
    :return: Tuple containing the average edge weights
    """
    # Initialize empty statistics
    statistics = [0.0, 0.0, 0.0]

    for bac in bacteria_pop:
        statistics[0] += bac.get_edge(f"{sugar}", f"transported_{sugar}")["weight"]
        statistics[1] += bac.get_edge(f"transported_{sugar}", f"enz_{sugar}_complex")["weight"]
        statistics[2] += bac.get_edge(f"enz_{sugar}_complex", "atp")["weight"]
    try:
        return tuple(map(lambda x: round(x[0] / len(bacteria_pop), 5), statistics))
    except TypeError as e:
        return tuple(map(lambda x: round(x / len(bacteria_pop), 5), statistics))

with open('records.tsv', 'w') as file:
    writer = csv.writer(file, delimiter="\t")

    for i in range(10000):
        try:
            # Progress simulator by 1 generation
            simulator.progress()

            # Obtain population statistics
            generation = i
            population_size = len(simulator.bacteria)
            food_availability = simulator.food_generator.food
            glucose_stats = sugar_stats(simulator.bacteria, "glucose")
            lactose_stats = sugar_stats(simulator.bacteria, "lactose")
            sucrose_stats = sugar_stats(simulator.bacteria, "sucrose")

            # Write stats to file
            writer.writerow([generation, population_size, *tuple(map(lambda pair: pair[1], sorted(list(food_availability.items()), key=lambda x: x[0]))),*glucose_stats, *lactose_stats, *sucrose_stats])
            file.flush()

            if i % 500 == 0 and i != 0:
                print("Generation: {}".format(i))
                print("Food available: {}".format(simulator.food_generator.food))
                print("Population size: {}".format(len(simulator.bacteria)))

                print("Glucose -> Transported Glucose {0:.5f} Transported_Glucose -> Enz_Glucose_Complex {1:.5f} Enz_Glucose_Complex -> ATP {2:.5f}".format(glucose_stats[0], glucose_stats[1], glucose_stats[2]))
                print("Lactose -> Transported Lactose {0:.5f} Transported_Lactose -> Enz_Lactose_Complex {1:.5f} Enz_Lactose_Complex -> ATP {2:.5f}".format(lactose_stats[0], lactose_stats[1], lactose_stats[2]))
                print("Sucrose -> Transported Sucrose {0:.5f} Transported_sucrose -> Enz_Sucrose_Complex {1:.5f} Enz_Sucrose_Complex -> ATP {2:.5f}".format(sucrose_stats[0], sucrose_stats[1], sucrose_stats[2]))

                if input("Continue? (enter \"n\" to terminate) ") == "n":
                    print(f"\nSimulation ended as generation {i}")
                    break

            time.sleep(0.1)

        except ZeroDivisionError as e:
            print(f"\nSimulation ended as generation {i}")
            break

