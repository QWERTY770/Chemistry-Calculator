from calculator import WaterSystem
from substances import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    x = []
    y = []
    system = WaterSystem()
    # sub1 = StrongAcid("HCl", SingleSpeciesGroup("HCl Species", Species("Cl-", -1)))
    sub1 = Acid.from_pKa("H3PO4", 2.12, 7.21, 12.67)
    sub2 = StrongBase("NaOH", SingleSpeciesGroup("NaOH Species", Species("Na+", 1)))
    system.add(Decimal("1.00"), sub1)
    for i in range(201):
        print(i, system.solve_pH())
        x.append(i/100)
        y.append(system.solve_pH())
        system.add(Decimal("0.01"), sub2)
    plt.title("1.00 mol/L H3PO4 + 1.00 mol/L NaOH * n")
    # plt.title("1.00 mol/L HCl + 1.00 mol/L NaOH * n")
    plt.xlabel("n(NaOH) added")
    plt.ylabel("pH")
    plt.plot(x, y, color="skyblue", linestyle="-", linewidth="2")
    plt.show()
