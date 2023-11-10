import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
import Constants

def processData(file_name):
    ''' data --> vectors '''
    file = open(file_name, 'r')
    processing = False
    x, y, data = [], [], []

    for line in file:
        if 'A004USERDEFINED' in line:
            processing = True
            continue
        
        if processing == True:
            cleaned = line.split()[1:]
            x = x + [int(cleaned[0]) + i for i in range(5)]
            y.extend(cleaned[1:])
        
    x = np.array(x)
    for i, item in enumerate(y):
        y[i] = int(item)
    y = np.array(y) # counts
    
    return x, y

def fit_gaussian(x, H, A, mu, sigma):
    ''' gaussian fit '''
    return A*np.exp((-(x-mu)**2)/(2*(sigma**2))) + H

def get_x_val(energy, m, b): 
    '''get the x index that matches the energy'''
    return int((energy-b)/(m))

def get_energy(channel):
    return (Constants.CALIB_SLOPE*channel) + Constants.CALIB_INT

def get_gaussian(start_energy, end_energy, x, y, m, b):
    start_range, end_range = get_x_val(start_energy, m, b), get_x_val(end_energy, m, b)
    x = x[start_range:end_range]
    y = y[start_range:end_range]
    mean = sum(x*y)/sum(y)              
    sigma = sum(y*(x-mean)**2)/sum(y) 
    print(min(y), max(y), (start_energy+end_energy)/2, sigma)
    return get_mu_sigma(curve_fit(fit_gaussian, x, y, p0=[min(y), max(y), mean, sigma]))

def get_mu_sigma(results):
    ''' get the mean and sigma values for the gaussian from the curve fit '''
    return results[0][2], results[0][3]

def get_fwhm(sigma):
    return 2 * math.sqrt(2*math.log(2)) * sigma

def plot_gaussian(start_energy, end_energy, x, y, params):
    ''' plot the gaussian '''
    fit_x = np.linspace(np.min(x_hist),np.max(x_hist),500)
    plt.plot(fit_x, fit_gaussian(fix_x, *params), 'r.:', label='gaussian fit')
    plt.legend()

def get_peak(y, l, r):
    l = get_x_val(l, Constants.CALIB_SLOPE, Constants.CALIB_INT)
    r = get_x_val(r, Constants.CALIB_SLOPE, Constants.CALIB_INT)
    print(l, r)

    count = 0
    sum = 0
    for i in range(l, r+1):
        energy = (Constants.CALIB_SLOPE * i) + Constants.CALIB_INT
        sum = sum + (y[i] * energy)
        count += y[i]
    
    return sum/count

def get_width_of_foil(observed, original, dE_e, dE_N):
    return (original - observed) / (dE_e + dE_N)

def get_error(observed, original, dE_e, dE_N):
    observed_e = observed*Constants.RESOLUTION
    original_e = original * Constants.RESOLUTION

    return round((((observed_e ** 2) + (original_e ** 2)) ** (0.5)) / (dE_e + dE_N), 3)

if __name__ == '__main__':
    path = '../data/'
    data_file = 'Al-Foil-10-26-2023.IEC'
    calib = '47-5V-10-26-2023.IEC'
    plot = True # plot option

    # get the and transform the vectors
    x, y = processData(path + data_file)
    xc, yc = processData(path + calib)

    m = Constants.CALIB_SLOPE
    b = Constants.CALIB_INT
    x = (m*x) + b # keV
    xc = (m*xc) + b
    
    # get resolution from 3182 keV peak
    GD_interval = [3165, 3193]
    AL_interval = [5460, 5489]

    r = get_gaussian(GD_interval[0], GD_interval[1], xc, yc, m, b)
    fwhm = get_fwhm(r[1])
    GD_148 = r[0]
    AM_241 = (Constants.AM_PEAK * m) + b
    print(f"FWHM: {fwhm}")
    print(f"Resolution: {(fwhm//m)/(get_x_val(r[0], m, b))}")
    print(f"Mean: {r[0]}, {AM_241}")

    # calculate widths
    # I found the peaks manually since I was unable to fit gaussians on GD and AM non-foil graphs
    HAV_23_GD = get_energy(471)
    HAV_23_AM = get_energy(1176)
    HAV_46_GD = get_energy(111)
    HAV_46_AM = get_energy(943)
    AL_GD = get_energy(299)
    AL_AM = get_energy(1079)
    NO_FOIL_GD = get_energy(808)
    NO_FOIL_AM = get_energy(1400)

    width = get_width_of_foil(AL_GD, NO_FOIL_GD, Constants.dE_ALUM_E_3180, Constants.dE_ALUM_N_3180)
    print(f"GD calculated width (AL): {round(width, 3)} +/- {get_error(AL_GD, NO_FOIL_GD, Constants.dE_ALUM_E_3180, Constants.dE_ALUM_N_3180)} microns")
    
    width = get_width_of_foil(AL_AM,  NO_FOIL_AM, Constants.dE_ALUM_E_5490, Constants.dE_ALUM_N_5490)
    print(f"AM calculated width (AL): {round(width, 3)} +/- {get_error(AL_AM,  NO_FOIL_AM, Constants.dE_ALUM_E_5490, Constants.dE_ALUM_N_5490)} microns")
    
    width = get_width_of_foil(HAV_23_GD, NO_FOIL_GD, Constants.dE_HAV_E_3180, Constants.dE_HAV_N_3180)
    print(f"GD calculated width (HAV 2.3): {round(width, 3)} +/- {get_error(HAV_23_GD, NO_FOIL_GD, Constants.dE_HAV_E_3180, Constants.dE_HAV_N_3180)} microns")

    width = get_width_of_foil(HAV_23_AM, NO_FOIL_AM, Constants.dE_HAV_E_5490, Constants.dE_HAV_N_5490)
    print(f"AM calculated width (HAV 2.3): {round(width, 3)} +/- {get_error(HAV_23_AM, NO_FOIL_AM, Constants.dE_HAV_E_5490, Constants.dE_HAV_N_5490)} microns")

    width = get_width_of_foil(HAV_46_GD, NO_FOIL_GD, Constants.dE_HAV_E_3180, Constants.dE_HAV_N_3180)
    print(f"GD calculated width (HAV 4.6): {round(width, 3)} +/- {get_error(HAV_46_GD, NO_FOIL_GD, Constants.dE_HAV_E_3180, Constants.dE_HAV_N_3180)} microns")

    width = get_width_of_foil(HAV_46_AM, NO_FOIL_AM, Constants.dE_HAV_E_5490, Constants.dE_HAV_N_5490)
    print(f"AM calculated width (HAV 4.6): {round(width, 3)} +/- {get_error(HAV_46_AM, NO_FOIL_AM, Constants.dE_HAV_E_5490, Constants.dE_HAV_N_5490)} microns")

    if plot:
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Count")
        ax.set_title("Calibration")
        plt.bar(xc,yc)
        plt.show()