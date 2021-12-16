# Standard Packages
from copy import deepcopy
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
from main.enums import SaveType
from main.variables import InputVariables, OutputVariables
from plots.styles import single as plotlayout
from plots.colors import mmm as plotcolors


class VarData:
    def __init__(self, xvar_name, yvar_name, label, runid, scan_num, var_to_scan=None, scan_factor=None, rho_value=None):
        self.xvar_name = xvar_name
        self.yvar_name = yvar_name
        self.label = label
        self.runid = runid
        self.scan_num = scan_num
        self.var_to_scan = var_to_scan
        self.scan_factor = scan_factor
        self.rho_value = rho_value
        self.xvar = None
        self.yvar = None


def plot_variable_data(data_list, title):
    plt.figure()

    for data in data_list:
        plt.plot(data.xvar.values, data.yvar.values, label=data.label)

    # Uses data from the last item in data_list
    plt.xlim(min([data.xvar.values.min() for data in data_list]),
             max([data.xvar.values.max() for data in data_list]))
    plt.xlabel(f'{data.xvar.label} {data.xvar.units_label}')
    plt.ylabel(f'{data.yvar.label} {data.yvar.units_label}')
    plt.legend()
    plt.title(title)
    plt.show()


def load_variable_data(data_list):

    plotlayout.init()
    plotcolors.init()

    for data in data_list:
        input_vars = InputVariables()
        output_vars = OutputVariables()

        args = (data.runid, data.scan_num, data.var_to_scan, data.scan_factor, data.rho_value)
        input_vars.load_from_csv(SaveType.INPUT, *args)
        input_vars.load_from_csv(SaveType.ADDITIONAL, *args)
        output_vars.load_from_csv(SaveType.OUTPUT, *args)

        var_list = [input_vars, output_vars]
        for v in var_list:
            v.set_rho_values()
            if hasattr(v, data.xvar_name):
                data.xvar = deepcopy(getattr(v, data.xvar_name))
            if hasattr(v, data.yvar_name):
                data.yvar = deepcopy(getattr(v, data.yvar_name))


if __name__ == '__main__':

    title, data_list = 'Collisionality', [
        VarData('rho', 'nuei', '120968A02 (High)', '120968A02', scan_num=1),
        VarData('rho', 'nuei', '120982A09 (Med.)', '120982A09', scan_num=1),
        VarData('rho', 'nuei', '129041A10 (Low)', '129041A10', scan_num=1),
    ]

    # title, data_list = 'Effective Charge', [
    #     VarData('rho', 'zeff', '120968A02 (High)', '120968A02', scan_num=1),
    #     VarData('rho', 'zeff', '120982A09 (Med.)', '120982A09', scan_num=1),
    #     VarData('rho', 'zeff', '129041A10 (Low)', '129041A10', scan_num=1),
    # ]

    # title, data_list = 'Collisionality', [
    #     VarData('rho', 'nuei', r'$\nu_\mathrm{ei}$', '120968A02', scan_num=1, var_to_scan='nuei', scan_factor=1),
    #     VarData('rho', 'nuei', r'$\nu_\mathrm{ei}\times 2$', '120968A02', scan_num=1, var_to_scan='nuei', scan_factor=2),
    #     VarData('rho', 'nuei', r'$\nu_\mathrm{ei}\times 3$', '120968A02', scan_num=1, var_to_scan='nuei', scan_factor=3),
    # ]

    # title, data_list = 'Growth Rate', [
    #     VarData('rho', 'gmaETGM', '120968A02 (High)', '120968A02', scan_num=1),
    #     VarData('rho', 'gmaETGM', '120982A09 (Med.)', '120982A09', scan_num=1),
    #     VarData('rho', 'gmaETGM', '129041A10 (Low)', '129041A10', scan_num=1),
    # ]

    title, data_list = r'Growth Rate $(\rho = 0.72)$', [
        VarData('nuei', 'gmaETGM', '120968A02 (High)', '120968A02', scan_num=1, var_to_scan='nuei', rho_value=0.72),
        VarData('nuei', 'gmaETGM', '120982A09 (Med.)', '120982A09', scan_num=1, var_to_scan='nuei', rho_value=0.72),
        VarData('nuei', 'gmaETGM', '129041A10 (Low)', '129041A10', scan_num=1, var_to_scan='nuei', rho_value=0.72),
    ]

    # title, data_list = 'Growth Rate (120968A02)', [
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $\nu_\mathrm{ei} \times 3$', '120968A02', scan_num=1, var_to_scan='nuei', scan_factor=3),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 3$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Growth Rate (120968A02, High)', [
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 0.5$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 1$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 2$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 3$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Growth Rate (120982A09, Med.)', [
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 0.5$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 1$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 2$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 3$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Growth Rate (129041A10, Low)', [
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 0.5$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 1$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 2$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'gmaETGM', r'$\gamma_\mathrm{etgm}$ for $Z_\mathrm{eff} \times 3$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Diffusivity (120968A02, High)', [
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 0.5$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 1$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 2$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 3$', '120968A02', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Diffusivity (120982A09, Med.)', [
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 0.5$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 1$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 2$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 3$', '120982A09', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    # title, data_list = 'Diffusivity (129041A10, Low)', [
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 0.5$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=0.5),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 1$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=1),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 2$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=2),
    #     VarData('rho', 'xteETGM', r'$\chi_\mathrm{e, etgm}$ for $Z_\mathrm{eff} \times 3$', '129041A10', scan_num=2, var_to_scan='zeff', scan_factor=3),
    # ]

    load_variable_data(data_list)
    plot_variable_data(data_list, title)
