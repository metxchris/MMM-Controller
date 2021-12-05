from enum import Enum

class PlotType(Enum):
    NONE = 0
    INPUT = 1
    OUTPUT = 2
    COMPARED = 3
    ADDITIONAL = 4

class ShotType(Enum):
    NONE = 0
    NSTX = 1
    DIII_D = 2
