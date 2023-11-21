from .conduction_velocity import ConductionVelocity
from .divergence import Divergence

class Analyse:

    def __init__(self, case):
        self.conduction_velocity = ConductionVelocity(case)
        self.divergence = Divergence(case)