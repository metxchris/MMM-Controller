from enum import Enum

class PlotType(Enum):
    INPUT = 1
    OUTPUT = 2
    COMPARED = 3
    ADDITIONAL = 4

class ShotType(Enum):
    NSTX = 1
    D3D = 2
