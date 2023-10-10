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
    
    return x, y

def fit_gaussian(x, H, A, mu, sigma):
    ''' gaussian fit '''
    return A*np.exp(-(x-mu)**2/(2*sigma**2)) + H

def get_x_val(energy): 
    '''get the x index that matches the energy'''
    return int((energy-1.36)/(0.243))

def get_gaussian(start_energy, end_energy, x, y):
    start_range, end_range = get_x_val(start_energy), get_x_val(end_energy)
    
    x = x[start_range:end_range]
    y = y[start_range:end_range]

    mean = sum(x*y)/sum(y)                  
    sigma = sum(y*(x-mean)**2)/sum(y) 

    return get_mu_sigma(curve_fit(fit_gaussian, x, y, p0=[min(y), max(y), mean, sigma]))

def get_mu_sigma(results):
    ''' get the mean and sigma values for the gaussian from the curve fit '''
    return results[0][2], np.diag(results[1])[2]

def plot_gaussian(start_energy, end_energy, x, y, params):
    ''' plot the gaussian '''
    fit_x = np.linspace(np.min(x_hist),np.max(x_hist),500)
    plt.plot(fit_x, fit_gaussian(fix_x, *params), 'r.:', label='gaussian fit')
    plt.legend()

def get_moment(factor):
    ''' calculate the moment of inertia '''
    h = 4.135*(10**-15)*(10**-3) #keV * s

    return((2*factor)/(h**2))

if __name__ == '__main__':
    data_file = 'Ho166_9-28-23.IEC'
    
    # get the and transform the vectors
    x, y = processData(data_file)
    x = np.array(x)
    for i, item in enumerate(y):
        y[i] = int(item)
    y = np.array(y) # counts
    x = (0.243*x) + 1.36 # keV

    # calculate gaussian curve fits for corresponding gamma rays
    intervals = [[361.8, 367.7], [276.6, 281.81], [180.9, 185.92], [78, 81.5]]
    intervals.reverse()
    results = []
    for i in intervals:
        r = get_gaussian(i[0], i[1], x, y)
        results.append([r[0], r[1]])
    
    # consolidate results into mean and standard deviation arrays
    levels, std = [0] * len(intervals), [0] * len(intervals)
    levels[0] = results[0][0]
    std[0] = results[0][1]
    for i in range(1, len(results)):
        levels[i] = levels[i-1] + results[i][0]
        std[i] = std[i-1] + results[i][1]

    # match levels to factors from red book and get moment of inertia
    levels,std = np.array(levels),np.array(std)
    factors = np.array([6, 20, 42, 72])
    print('moment: {0} +- {1}'.format(get_moment(factors), get_moment(std)))
    print('level/factors: {0}'.format(levels/factors))
    
    ax = plt.gca()
    ax.set_xlim(min(x), 2000)
    ax.set_ylim(min(y), 50000)
    plt.bar(x,y)
    plt.show()
