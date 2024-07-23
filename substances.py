from species import *
from typing import Iterable

zero = Decimal("0")
inf = Decimal("Infinity")
ninf = Decimal("-Infinity")
neutral = Decimal("1E-7")
Kw = Decimal("1E-14")


def pKa(f) -> Decimal:
    return Decimal("10") ** -Decimal(str(f))


class Substance:
    def __init__(self, name: str, species: Iterable[tuple[SpeciesGroup, int]] = None):
        self.name = name
        self.species = species if species is not None else []  # type: Iterable[tuple[SpeciesGroup, int]]


class SingleSubstance(Substance):
    def __init__(self, name: str, species: SingleSpeciesGroup, num: int = 1):
        super().__init__(name, ((species, num),))


class Acid(Substance):
    def __init__(self, name: str, species: AcidSpeciesGroup, num: int = 1):
        super().__init__(name, ((species, num),))

    @classmethod
    def from_Ka(cls, name: str, *Kas: Decimal):
        sp = []
        for i in range(len(Kas)+1):
            sp.append(Species(f"{name}_species_{i}", -i))
        return Acid(name, AcidSpeciesGroup(f"{name} Species", sp, *Kas))

    @classmethod
    def from_pKa(cls, name: str, *pKas):
        Kas = tuple(map(pKa, pKas))
        return Acid.from_Ka(name, *Kas)


StrongAcid = SingleSubstance
StrongBase = SingleSubstance
