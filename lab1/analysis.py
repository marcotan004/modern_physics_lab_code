import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit

def processData(file_name):
    file = open(file_name, 'r')
    processing = False
    x, y, data = [], [], []

    for line in file:
        if 'A004USERDEFINED' in line:
            processing = True
            continue
        
        if processing == True:
            cleaned = ' '.join(line.split()).replace('A004 ', '').split()
            x = x + [float(cleaned[0])]
            vals = [float(i) for i in cleaned[1:]]

            y = y + [sum(vals)]
    
    return x, y

def getTimeData(x, y, n):
    # was used to find peaks in time_data
    maxes = [(0,0)] * n 

    for i, num in enumerate(y):
        if num >= maxes[0][0]:
            maxes[0] = (num, x[i])
            maxes.sort(key=lambda a:a[0])
    
    return maxes

def getHalfLife(b):
    return math.log(2)/b


if __name__ == '__main__':
    data_file = '9-21-2023.IEC'
    
    time_file = '09-05-2023-Time-Calibration.IEC'
    start = 228
    end = 541
    
    # 330 channels --> 0.16 microseconds
    x, y = processData(data_file)
    
    # time calibration factor found by looking at the graph from time_data
    channelTimeFactor = 0.16/330.0 #* 10**-6
    times = [round(num*channelTimeFactor,5) for num in x]
    
    t0 = times[start]

    for i, num in enumerate(y):
        print(i, num)
    print(t0)

    # exponential function
    def func(x, a, b, c):
        return a * np.exp(-b * (x-t0)) + c

    popt, pcov = curve_fit(func, times[start:end], y[start:end])
    print('y = {0} * exp(-{1} * x) + {2}'.format(popt[0], popt[1], popt[2]))

    # standard deviations of ['a', 'b', 'c'] in perr
    perr = np.diag(pcov)
    print(pcov)
    print('b = {0} += {1}'.format(popt[1], perr[1]))
    print('Halflife: {0} microseconds'.format(round(getHalfLife(popt[1]),4)))
    print(times[-1])
    plt.figure(1, figsize=(5,5))
    plt.suptitle('Experiment 2 Count Data')

    # set mins and maxes on graph
    ax = plt.gca()
    ax.set_xlim(min(times), max(times))
    ax.set_ylim(min(y), max(y) + 100)

    # to calculate the fit
    linsp = np.array(times)
    plt.scatter(times, y, s=1, label='data')
    plt.plot(linsp, func(linsp, *popt), label=r'${y = 505.6e^{-7.89(x-0.55)} + 0.72}$', c='red')
    #r'${y = 96.9e^{-7.47(x-0.61)} + 0.27}$'
    #label=f'y = {round(popt[0],1)} * exp(-{round(popt[1],2)} * x) + {round(popt[2],2)}'
    plt.xlabel(u'Time (\u03bcs)')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

    