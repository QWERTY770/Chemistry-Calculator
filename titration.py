from calculator import WaterSystem
from substances import *
import matplotlib.pyplot as plt


class Titration:
    """TODO"""
    def __init__(self, substances: Iterable[tuple[Substance, Decimal]], system, *args):
        self.system = system
        self.substances = substances

    def get_data(self):
        pass  # Not implemented


class AcidBaseTitration(Titration):
    """TODO"""
    def __init__(self, substances: Iterable[tuple[Substance, Decimal]], system, *args):
        super().__init__(substances, system, args)

    def get_data(self):
        pass


if __name__ == "__main__":
    x = []
    y = []
    system = WaterSystem()

    Na = SingleSpeciesGroup("Na+", Species("Na+", 1))
    EDTA_H6 = Species("H6-EDTA 2+", 2)
    EDTA_H5 = Species("H5-EDTA +", 1)
    EDTA_H4 = Species("H4-EDTA", 0)
    EDTA_H3 = Species("H3-EDTA -", -1)
    EDTA_H2 = Species("H2-EDTA 2-", -2)
    EDTA_H1 = Species("H-EDTA 3-", -3)
    EDTA_H0 = Species("EDTA 4-", -4)
    EDTA = AcidSpeciesGroup("EDTA Species", [EDTA_H6, EDTA_H5, EDTA_H4, EDTA_H3, EDTA_H2, EDTA_H1, EDTA_H0],
                            pKa(1.15), pKa(1.15), pKa(2.12), pKa(2.57), pKa(6.16), pKa(10.26))
    # M.A. Marini, W.J. Evans, R.L. Berger, Use of the twin-cell differential titration calorimeter for binding studies.
    # I. EDTA and its calcium complex, J. Biochem. Biophys. Methods 10 1985 273.

    sub1 = Substance("Na4EDTA", ((Na, 4), (EDTA, 1)))
    sub2 = StrongAcid("HCl", SingleSpeciesGroup("HCl Species", Species("Cl-", -1)))
    system.substances[sub2] = Decimal("0")
    system.add(Decimal("0.10"), sub1)
    v = Decimal("50")
    per = Decimal("1")
    for i in range(300):
        pH = system.solve_pH()
        print(i*per, pH, v, system.substances[sub1], system.substances[sub2])
        system.add(Decimal("0.10") * per / v, sub2)
        x.append(i*per)
        y.append(pH)
        system.set_volume(v / (v + per))
        v += per
    plt.title("0.10 mol/L Na4EDTA + 0.10 mol/L HCl")
    plt.xlabel("V(HCl) added")
    plt.ylabel("pH")
    plt.plot(x, y, color="skyblue", linestyle="-", linewidth="2")
    plt.show()
