from typing import List

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_ID = {aa: i for i, aa in enumerate(AA)}
ID_TO_AA = {i: aa for aa, i in AA_TO_ID.items()}

class Individual:

    counter = 0

    def __init__(self, sequence: str):
        self.id = str(Individual.counter)
        Individual.counter += 1

        self.X = self._aa2num(sequence)
        self.F: np.ndarray = None      # espacio de minimización (usado por MOEA/D)
        self.fitness: List[str] = None
        self.F_raw: np.ndarray = None  # valores originales sin transformar (para guardar)
        self.father: "Individual" = None
        self.pdb: str = None
        self.nmut: int = 0

    def _aa2num(self, sequence: str):
        return np.array([AA_TO_ID[aa] for aa in sequence], dtype=np.int32)

    def num2aa(self):
        return [ID_TO_AA[x] for x in self.X]

    def sequence(self):
        return "".join(self.num2aa())

    def copy(self):
        new = Individual(self.sequence())
        new.X = self.X.copy()
        new.F     = None if self.F     is None else self.F.copy()
        new.F_raw = None if self.F_raw is None else self.F_raw.copy()
        new.father = self.father
        new.pdb = self.pdb
        new.nmut = self.nmut
        return new

    def __repr__(self):
        return f"Individual( id={self.id}, seq='{self.sequence()}')"