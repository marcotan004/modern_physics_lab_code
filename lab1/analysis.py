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
    data_file = '09-19-2023.IEC'
    time_file ='09-05-2023-Time-Calibration.IEC'
    start = 228
    end = 500
    
    # 330 channels --> 0.16 microseconds
    x, y = processData(data_file)

    # time calibration factor found by looking at the graph from time_data
    channelTimeFactor = 0.16/330.0 #* 10**-6
    times = [round(num*channelTimeFactor,5) for num in x]
    
    t0 = times[start]

    # exponential function
    def func(x, a, b, c):
        return a * np.exp(-b * (x-t0)) + c

    popt, pcov = curve_fit(func, times[start:end], y[start:end])
    print('y = {0} * exp(-{1} * x) + {2}'.format(popt[0], popt[1], popt[2]))

    # standard deviations of ['a', 'b', 'c'] in perr
    perr = np.diag(pcov)
    print('b = {0} += {1}'.format(popt[1], perr[1]))
    print('Halflife: {0} microseconds'.format(round(getHalfLife(popt[1]),4)))

    # for checking if it is an exponential curve
    log_times = [math.log(i) for i in y[start:end]]

    plt.figure(1, figsize=(5,5))
    plt.suptitle('Count Data')

    # set mins and maxes on graph
    ax = plt.gca()
    ax.set_xlim(min(times), max(times))
    ax.set_ylim(min(y), max(y) + 100)

    # to calculate the fit
    linsp = np.array(times)
    plt.plot(linsp, func(linsp, *popt), times, y, label='fit')

    plt.xlabel('MicroSeconds')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

    