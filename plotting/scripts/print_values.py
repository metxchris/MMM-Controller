# Standard Packages
import sys; sys.path.insert(0, '../'), sys.path.insert(0, '../../')
import copy
import os

# 3rd Party Packages
import matplotlib.pyplot as plt
import numpy as np

# Local Packages
import settings
import modules.controls
import modules.constants
import modules.options
import modules.calculations as calculations
import modules.adjustments as adjustments
import modules.datahelper as datahelper
import modules.mmm as mmm
import modules.reshaper as reshaper
import modules.utils as utils
import modules.variables as variables
import modules.controls as controls
import plotting.modules.profiles as profiles
from modules.enums import ShotType, ScanType, ProfileType, SaveType


runid, input_time, scan_num = '129016A03', 0.46, 30  # wexb 0.5
# runid, input_time, scan_num = '129016A03', 0.46, 26  # wexb on
runid, input_time, scan_num = '129016A03', 0.46, 27  # wexb off
runid, input_time, scan_num = '129016A03', 0.46, 40  # wexb off, mi = mh
runid, input_time, scan_num = '129016A03', 0.46, 41  #0.5 fsa
# runid, input_time, scan_num = '129016A04', 0.49, 19

def print_value(name, vars):
    print(name + ',', f'{get_value(name, vars):.3e}')

def get_value(name, vars):
    return getattr(vars, name).values[r, options.time_idx][0]


options = modules.options.Options(runid=runid, scan_num=scan_num)
# input_vars, cdf_vars, __ = datahelper.initialize_variables(options)

input_vars = variables.InputVariables(options)
output_vars = variables.OutputVariables(options)
controls = controls.InputControls(options)
controls.load_from_csv()

input_vars.load_from_csv(SaveType.INPUT)
input_vars.load_from_csv(SaveType.ADDITIONAL)
output_vars.load_from_csv(SaveType.OUTPUT)

zckb = modules.constants.ZCKB

for rmina in [0.6, 0.7]:

    r = np.absolute(input_vars.rmina.values[:, options.time_idx] - rmina).argmin()

    print(f'Variable values for r/a = {input_vars.rmina.values[r, options.time_idx][0]:.3f}:')
    print('-' * 50)

    amin = input_vars.rmin.values[-1, options.time_idx][0]

    print_value('rho', input_vars)
    print_value('rmina', input_vars)
    print_value('rmin', input_vars)
    print(f'amin, {amin:.3e}')
    print_value('rmaj', input_vars)
    print_value('elong', input_vars)
    # print_value('gxi', input_vars)
    print_value('ne', input_vars)
    print_value('nf', input_vars)
    print_value('nh', input_vars)
    print_value('nz', input_vars)
    print_value('te', input_vars)
    print_value('ti', input_vars)
    print_value('btor', input_vars)
    print_value('bu', input_vars)
    print_value('q', input_vars)
    print_value('zeff', input_vars)
    print_value('zz', input_vars)
    # print_value('ai', input_vars)
    print_value('ah', input_vars)
    print_value('az', input_vars)
    print_value('wexb', input_vars)
    print('wexb (cs / a), 0' )
    print_value('gne', input_vars)
    print_value('gni', input_vars)
    print_value('gnh', input_vars)
    print_value('gnz', input_vars)
    rmaj = get_value('rmaj', input_vars)
    a_lte = amin / rmaj * get_value('gte', input_vars)
    print(f'a_lte, {a_lte:.3e}')
    print_value('gte', input_vars)
    print_value('gti', input_vars)
    print_value('gq', input_vars)
    print_value('gbu', input_vars)

    print_value('rhosu', input_vars)
    print_value('gyrfiu', input_vars)
    print_value('csound_a', input_vars)
    dv = 46 if rmina == 0.6 else 55
    print(f'V\', {dv:.3e}')
    # print_value('gaveETGM', output_vars)

    print_value('kyrhosETGM', output_vars)
    print('gma - wexb, 0')
    print_value('gmaETGM', output_vars)
    print_value('omgETGM', output_vars)
    print('gma - wexb (cs / a), 0')
    cs_a = get_value('csound_a', input_vars)
    gma = get_value('gmaETGM', output_vars)
    gman = gma / cs_a
    print(f'gman, {gman:.3e}')
    
    omg = get_value('omgETGM', output_vars)
    omgn = omg / cs_a
    print(f'omgn, {omgn:.3e}')
    print_value('xteETGM', output_vars)
    # print_value('xte2ETGM', output_vars)
    print_value('fte', output_vars)
    ftene = get_value('fte', output_vars) * get_value('ne', input_vars) * zckb
    print(f'Q_e, {ftene:.3e}')
    power = ftene * dv * amin**2 
    print(f'power, {power:.3e}')
    print()
    

