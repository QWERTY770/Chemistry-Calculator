from substances import *
from math import log10
from mpmath import mp, power
mp.prec = 53
mp.dps = 24
mp.eps = power(10, -30)


def check_solution(equation, mapping, solution, cH) -> None:
    print(equation)
    print(solution)
    for i in mapping:
        print(f"{i}={mapping[i].subs(cH, solution)}")


class WaterSystem:
    def __init__(self, substances: dict[Substance, Decimal] = None, Kw: Decimal = Kw):
        self.substances = substances if substances is not None else {}
        self.cH = Species("H+", 1)
        self.cOH = Species("OH-", -1)
        self.groups = {}  # type: dict[SpeciesGroup, Decimal]
        self.Kw = Kw
        for sub in self.substances:
            c = substances[sub]
            for spec in sub.species:
                self.groups[spec[0]] = spec[1] * c

    def add(self, c, sub: Substance):
        if isinstance(c, float):
            c = Decimal(c)
        if sub in self.substances:
            self.substances[sub] += c
        else:
            self.substances[sub] = c
        for spec in sub.species:
            if spec[0] in self.groups:
                self.groups[spec[0]] += spec[1] * c
            else:
                self.groups[spec[0]] = spec[1] * c

    def _solve(self, equation, initial_value=Decimal("1E-7")):
        tol = -30
        solution = nsolve(equation, self.cH, initial_value, maxsteps=100, tol=power(10, tol))
        if solution < 0:
            solution = nsolve(equation.subs(self.cH, self.Kw / self.cOH),
                              self.cOH, initial_value, maxsteps=100, tol=power(10, tol))
            if solution < 0:
                solution = nsolve(equation, (Decimal("1E-20"), 1),
                                  solver="bisect", verify=False, maxsteps=200, tol=power(10, tol))
                if solution < 0:
                    raise RuntimeError(f"Solve failed, equation={equation}, solution={solution}")
                return solution
            solution = self.Kw / solution
        return solution

    def solve(self, initial_value=Decimal("1E-7")):
        mapping = {}
        for spec in self.groups:
            c = self.groups[spec]
            mapping.update(spec.get_equations(c, self.cH))
        mapping[self.cOH] = self.Kw / self.cH
        charge_eq = self.cH - self.Kw / self.cH
        for i in self.groups:
            for j in i.species:
                charge_eq += mapping[j]*j.charge
        solution = self._solve(charge_eq, initial_value)
        # check_solution(charge_eq, mapping, solution, self.cH)
        return solution

    def solve_pH(self):
        return -log10(self.solve())


if __name__ == '__main__':
    system = WaterSystem()
    sub1 = Acid.from_pKa("H3PO4", 2.12, 7.21, 12.67)
    sub2 = Acid.from_pKa("CH3CHOHCOOH", 4.76)
    sub3 = StrongAcid("HCl", SingleSpeciesGroup("HCl Species", Species("Cl-", -1)))
    sub4 = StrongBase("NaOH", SingleSpeciesGroup("NaOH Species", Species("Na+", 1)))
    system.add(Decimal("1"), sub2)
    print(system.solve_pH())
