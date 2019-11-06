# importing libraries
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Color palette
# https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7 for colorpalette

fig = plt.figure()
plt.axis('off')
plt.title('Bacteria Simulator Live Monitoring', y=1.08)

# creating a subplot
pop_ax = fig.add_subplot(2, 2, 1)
glucose_ax = fig.add_subplot(2, 2, 2)
lactose_ax = fig.add_subplot(2, 2, 3)
sucrose_ax = fig.add_subplot(2, 2, 4)


def line_to_dict(line, delim="\t"):
    data = line.split(delim)
    return {
        "generation"    : float(data[0]),
        "population"    : float(data[1]),
        "glucose_amt"   : float(data[2]),
        "lactose_amt"   : float(data[3]),
        "sucrose_amt"   : float(data[4]),
        "glucose_trans" : float(data[5]),
        "glucose_enz"   : float(data[6]),
        "glucose_atp"   : float(data[7]),
        "lactose_trans" : float(data[8]),
        "lactose_enz"   : float(data[9]),
        "lactose_atp"   : float(data[10]),
        "sucrose_trans" : float(data[11]),
        "sucrose_enz"   : float(data[12]),
        "sucrose_atp"   : float(data[13])
    }

def animate(i):
    data = open('records.tsv', 'r').read()
    lines = data.split('\n')
    generation = []
    population = []
    glucose_amt = []
    lactose_amt = []
    sucrose_amt = []
    glucose_trans = []
    glucose_enz = []
    glucose_atp = []
    lactose_trans = []
    lactose_enz = []
    lactose_atp = []
    sucrose_trans = []
    sucrose_enz = []
    sucrose_atp = []

    for line in lines[-50:-1]:
        statistics = line_to_dict(line)  # Delimiter is tab
        generation.append(statistics["generation"])
        population.append(statistics["population"])
        glucose_amt.append(statistics["glucose_amt"])
        lactose_amt.append(statistics["lactose_amt"])
        sucrose_amt.append(statistics["sucrose_amt"])
        glucose_trans.append(statistics["glucose_trans"])
        glucose_enz.append(statistics["glucose_enz"])
        glucose_atp.append(statistics["glucose_atp"])
        lactose_trans.append(statistics["lactose_trans"])
        lactose_enz.append(statistics["lactose_enz"])
        lactose_atp.append(statistics["lactose_atp"])
        sucrose_trans.append(statistics["sucrose_trans"])
        sucrose_enz.append(statistics["sucrose_enz"])
        sucrose_atp.append(statistics["sucrose_atp"])

    # Figure settings
    pop_ax.clear()
    pop_ax.set_ylim([0, max([*population, *glucose_amt, *lactose_amt, *sucrose_amt]) * 1.1])
    pop_ax.set_title("Population & Food")

    glucose_ax.clear()
    glucose_ax.set_ylim([0, min(max([*glucose_trans, *glucose_enz, *glucose_atp]) * 1.1, 1)])
    glucose_ax.set_title("Glucose Pathway")

    lactose_ax.clear()
    lactose_ax.set_ylim([0, 1])
    lactose_ax.set_ylim([0, min(max([*lactose_trans, *lactose_enz, *lactose_atp]) * 1.1, 1)])
    lactose_ax.set_title("Lactose Pathway")

    sucrose_ax.clear()
    sucrose_ax.set_ylim([0, 1])
    sucrose_ax.set_ylim([0, min(max([*sucrose_trans, *sucrose_enz, *sucrose_atp]) * 1.1, 1)])
    sucrose_ax.set_title("Sucrose Pathway")

    pop_ax.plot(generation, population, "#D55E00", linestyle='dashed', label="Population")
    pop_ax.plot(generation, glucose_amt, "#0072B2", label="Glucse amt")
    pop_ax.plot(generation, lactose_amt, "#F0E442", label="Lactose amt")
    pop_ax.plot(generation, sucrose_amt, "#009E73", label="Sucrose amt")
    pop_ax.legend(loc="bottom left")

    glucose_ax.plot(generation, glucose_trans, "#D55E00", label="Transporter")
    glucose_ax.plot(generation, glucose_enz, "#0072B2", label="to Intermediate")
    glucose_ax.plot(generation, glucose_atp, "#009E73", label="to ATP")
    glucose_ax.legend(loc="right")

    lactose_ax.plot(generation, lactose_trans, "#D55E00", label="Transporter")
    lactose_ax.plot(generation, lactose_enz, "#0072B2", label="to Intermediate")
    lactose_ax.plot(generation, lactose_atp, "#009E73", label="to ATP")
    lactose_ax.legend(loc="right")

    sucrose_ax.plot(generation, sucrose_trans, "#D55E00", label="Transporter")
    sucrose_ax.plot(generation, sucrose_enz, "#0072B2", label="to Intermediate")
    sucrose_ax.plot(generation, sucrose_atp, "#009E73", label="to ATP")
    sucrose_ax.legend(loc="right")

ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()
