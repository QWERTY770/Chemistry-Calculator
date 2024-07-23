from decimal import *
from sympy import *


class Species(Symbol):
    def __new__(cls, formula, *args):
        return super().__new__(cls, formula, real=True)

    def __init__(self, formula: str, charge: int = 0):
        self.formula = formula
        self.charge = charge

    def __hash__(self):
        return hash((self.formula, self.charge))

    def as_charge(self):
        return self.charge * self


class SpeciesGroup:
    def __init__(self, name: str, species: list[Species] = None):
        self.name = name
        self.species = species

    def __hash__(self):
        return hash((self.name, tuple(self.species)))

    def get_equations(self, c, cH: Species) -> dict:
        pass


class SingleSpeciesGroup(SpeciesGroup):
    def __init__(self, name: str, species: Species):
        if not isinstance(species, Species):
            raise RuntimeWarning(f"Invalid single species group, species={species}!")
        super().__init__(name, [species])

    def get_equations(self, c, cH: Species) -> dict:
        c = Float(c)
        return {self.species[0]: c}


class AcidSpeciesGroup(SpeciesGroup):
    def __init__(self, name: str, species: list[Species] = None, *Ka: Decimal):
        self.Ka = Ka
        self.protons = len(Ka)
        super().__init__(name, species)

    def get_equations(self, c, cH: Species) -> dict:
        if self.protons <= 0:
            raise RuntimeError(f"Acid {self.name} is not a valid acid!")
        equations = {}
        terms = [cH ** self.protons]
        Ka_part = 1
        for i, j in enumerate(self.Ka):
            Ka_part *= j
            terms.append(cH ** (self.protons - i - 1) * Ka_part)
        deno = sum(terms)  # denominator
        for i, j in enumerate(terms):
            equations[self.species[i]] = c * j / deno
        return equations
