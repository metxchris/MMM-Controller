# Standard Packages
from enum import Enum


class PlotType(Enum):
    NONE = 0
    INPUT = 1
    OUTPUT = 2
    COMPARED = 3
    ADDITIONAL = 4

class DataType(Enum):
    NONE = 0
    INPUT = 1
    OUTPUT = 2
    CONTROL = 3

class ScanType(Enum):
    NONE = 0
    VARIABLE = 1
    CONTROL = 2

class ShotType(Enum):
    NONE = 0
    NSTX = 1
    DIII_D = 2
