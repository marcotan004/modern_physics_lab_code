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

def get_moment(energy, factor):
    ''' calculate the moment of inertia '''
    h = 6.582*(10**-16) #eV * s
    energy = energy*1000# keV -> eV

    return (h**2/(2*energy))*factor # (ev * s)**2 / (eV) = eV * (s ** 2) 

def get_moment_uncertainty(std, energy, factor):
    h = 6.582*(10**-16) #eV * s
    return (factor)*(h**2)*(1/2)*(std/energy)

def expected_energy(moment, factor):
    h = 6.582*(10**-16) 

    return ((h**2)/(2*(moment))) * factor

def phi_rigid(M, R, phi):
    ''' calculate beta of the rigid body model '''
    # phi in eV * (s ** 2) 
    # convert to kg * m^2
    conv = 1.602 * (10**-19)
    phi = phi * conv

    # calculate beta
    numerator = 5*phi
    denominator = 0.31 * 2 * M * R

    return (numerator/denominator) - 1 # returns beta

def phi_fluid(M, R, phi):
    # phi in eV * (s ** 2) 
    # convert to kg * m^2
    conv = 1.602 * (10**-19)
    phi = phi * conv

    # calculate beta
    numerator = 8 * phi * math.pi 
    denominator = 9 * M * (R ** 2)
    return math.sqrt(numerator/denominator)

if __name__ == '__main__':
    data_file = 'Ho166_9-28-23.IEC'
    plot = False # plot option

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
    moments = get_moment(levels, factors)
    uncert = get_moment_uncertainty(std, levels, factors)

    print('moment of inertia (eV * s**2): {0} +- {1}'.format(moments, uncert))
    print('mean moment of inertia: {0}'.format(np.mean(moments)))
    print('level/factors: {0}'.format(levels/factors))
    
    # calculate beta values using mean phi
    M = 165.93228 * (1.6605 * 10**-27) # atomic mass of Ho-166 (kg)
    R = (1.2 *(10 ** -15)) * (166 ** (1/3))
    phi = np.mean(moments)
    beta_rigid = round(phi_rigid(M, R, phi), 2)
    beta_fluid = round(phi_fluid(M, R, phi), 2)
    print('Beta of rigid body model: {0} | Beta of fluid model: {1}'.format(beta_rigid, beta_fluid))

    if plot:
        ax = plt.gca()
        ax.set_xlim(min(x), 2000)
        ax.set_ylim(min(y), 50000)
        plt.bar(x,y)
        plt.show()
