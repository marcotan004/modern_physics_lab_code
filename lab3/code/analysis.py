import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit

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

def get_gaussian(start_energy, end_energy, x, y, m, b):
    start_range, end_range = get_x_val(start_energy, m, b), get_x_val(end_energy, m, b)
    x = x[start_range:end_range]
    y = y[start_range:end_range]
    mean = sum(x*y)/sum(y)              
    sigma = sum(y*(x-mean)**2)/sum(y) 

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

if __name__ == '__main__':
    path = '../data/'
    data_file = '23-Foil-10-31-2023.IEC'
    calib = '47-5V-10-26-2023.IEC'
    plot = True # plot option

    # get the and transform the vectors
    x, y = processData(path + data_file)
    xc, yc = processData(path + calib)

    m = 3.883
    b = 37.4
    x = (m*x) + b # keV
    xc = (m*xc) + b
    
    # get resolution from 3182 keV peak
    interval = [3165, 3193]
    r = get_gaussian(interval[0], interval[1], xc, yc, m, b)
    fwhm = get_fwhm(r[1])
    print(f"FWHM: {fwhm}")
    print(f"Resolution: {(fwhm//m)/(get_x_val(r[0], m, b))}")

    if plot:
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Count")
        ax.set_title("Calibration")
        plt.bar(xc,yc)
        plt.show()