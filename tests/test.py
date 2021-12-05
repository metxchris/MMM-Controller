# Standard Packages
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np 
import matplotlib.pyplot as plt
import scipy.ndimage
# Local Packages
from main.options import Options

def movingaverage_test():
    window_width = 3
    VPOL = [-499.268,-499.268,-1160.13,-1906.77,-2639.36,-3271.14,-3839.17,-4327.92,-4721.99,-5020.67,-5235.96,-5377.69,-5451.62,-5465.43,-5424.32,-5330.84,-5188.84,-5004.10,-4780.58,-4519.62,-4230.65,-3927.62,-3599.96,-3244.53,-2899.74,-2584.62,-2283.50,-1992.06,-1723.03,-1486.91,-1282.03,-1098.73,-937.240,-815.939,-743.808,-712.363,-713.436,-740.849,-788.835,-856.872,-948.787,-1054.12,-1142.63,-1195.56,-1220.26,-1239.22,-1362.68,-1724.38,-2195.38,-2018.60,-1446.27]
    cumsum_vec = np.cumsum(np.insert(VPOL, 0, 0)) 
    VPOL_smoothed = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    x = range(len(VPOL))
    x2 = range(1, VPOL_smoothed.size + 1)
    plt.figure()
    plt.plot(x, VPOL, x2, VPOL_smoothed)
    plt.show()

def smoothing_test():
    VPOL = [-499.268,-499.268,-1160.13,-1906.77,-2639.36,-3271.14,-3839.17,-4327.92,-4721.99,-5020.67,-5235.96,-5377.69,-5451.62,-5465.43,-5424.32,-5330.84,-5188.84,-5004.10,-4780.58,-4519.62,-4230.65,-3927.62,-3599.96,-3244.53,-2899.74,-2584.62,-2283.50,-1992.06,-1723.03,-1486.91,-1282.03,-1098.73,-937.240,-815.939,-743.808,-712.363,-713.436,-740.849,-788.835,-856.872,-948.787,-1054.12,-1142.63,-1195.56,-1220.26,-1239.22,-1362.68,-1724.38,-2195.38,-2018.60,-1446.27]
    x = range(len(VPOL))
    fig, ax = plt.subplots()
    ax.plot(x, VPOL, 'b')
    ax.plot(x, scipy.ndimage.gaussian_filter(VPOL, 1), 'r')
    plt.show()

def plot_gradient():
    gTE = [0,0.893,0.976,1.091,1.276,1.470,1.694,1.921,2.149,2.374,2.593,2.803,3.007,3.205,3.399,3.587,3.770,3.947,4.117,4.282,4.442,4.596,4.748,4.897,5.045,5.189,5.326,5.460,5.587,5.712,5.826,5.932,6.034,6.122,6.211,6.304,6.422,6.623,6.998,7.696,8.962,11.07,14.39,19.24,26.31,36.71,49.82,64.98,86.95,100,100]
    x = range(len(gTE))
    fig, ax = plt.subplots()
    ax.plot(x, gTE, 'b')
    ax.plot(x, scipy.ndimage.gaussian_filter(gTE, 1), 'r')
    plt.show()

def plot_etai(vars):
    ETAI = [None, 0.04789,1.10650,1.18148,1.29022,1.38803,1.48440,1.56889,1.64440,1.71295,1.77916,1.84783,1.92343,2.00964,2.10964,2.22678,2.36572,2.53419,2.74499,3.01798,3.38332,3.88659,4.59630,5.61909,7.13489,9.48201,13.4705,21.9128,57.0906,-77.035,-20.586,-11.047,-7.2936,-5.4679,-4.5580,-4.2138,-4.3440,-5.1090,-7.4521,-21.304,15.3831,4.98349,2.81117,1.90058,1.43372,1.14140,0.96853,0.82807,0.73724,0.65212,0.79173]
    gTI = [0,1.320,1.455,1.639,1.933,2.240,2.593,2.949,3.312,3.669,4.018,4.355,4.678,4.987,5.280,5.550,5.795,6.007,6.180,6.312,6.395,6.426,6.404,6.324,6.187,5.991,5.739,5.446,5.125,4.800,4.486,4.209,3.997,3.862,3.835,3.933,4.177,4.600,5.230,6.110,7.313,8.913,11.01,13.71,17.27,22.24,27.57,32.58,37.41,46.62,79.17]
    gNI = [0,1.2450,1.3599,1.5115,1.7502,1.9856,2.2443,2.4862,2.7114,2.9096,3.0781,3.2132,3.3151,3.3842,3.4202,3.4201,3.3822,3.3035,3.1804,3.0150,2.8090,2.5693,2.3092,2.0431,1.7859,1.5473,1.3323,1.1389,0.9578,0.7787,0.5917,0.3938,0.1897,-0.009,-0.189,-0.328,-0.408,-0.399,-0.254,0.1068,0.8176,2.0819,4.1656,7.4225,12.219,19.620,28.576,39.414,50.799,71.542,100]
    print(vars)
    rho = vars.rho.values[:,-1]
    plt.figure()
    # plt.plot(rho, vars.gti.values[:,0], rho, gTI)
    # plt.plot(rho, vars.gni.values[:,0], rho, gNI)
    plt.plot(rho, vars.etai.values[:,-1], rho, ETAI)
    plt.show()

def test_options():
    print(Options.runid, Options.runid)
    # Options.time=7

if __name__ == '__main__':
    ...


    '''
    REGEX Search Testing
    '''
    # @units_label.setter
    # def units_label(self):
    #     print('units_label')
    #     # Convert self.units into LaTeX format
    #     search_strs = ['\^\d', '\^\-\d', '\*\*\d', '\*\*\-\d']
    #     for s in search_strs:
    #         search_result = re.compile(s).search(self._units)
    #         print(self.__name__, search_result)

    # search_strs = ['\^\d', '\^\-\d', '\*\*\d', '\*\*\-\d']
    # for s in search_strs:
    #     search_result = re.compile(s).search(self._units)
    #     if search_result is not None:
    #         search_str = search_result.group()
    #         number_in_str = re.compile('\d').search(search_str).group()
    #         print(number_in_str)
    #         print(self._units.replace(search_result.group(), r'\$\^\{' + number_in_str + r'\}\$'))
    #     print(self.name, search_result)
