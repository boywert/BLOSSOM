import numpy
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

from math import *
class parameters(object):
    z = 8.064
    binfields = 6
    max_resolve = 100.
    config = {}
    
def readconfig(configfile):
    f = open(configfile)
    line = f.read().splitlines()
    config = {}
    for i in range(len(line)):
        column = line[i].split()
        if len(column) > 1:
            if column[0][0] != '#':
                if column[0][-5:]=='_path' or column[0][-5:]=='_file':
                    config[column[0]] = column[1]
                    parameters.config[column[0]] = column[1]
                else:
                    if column[1].find('.') > -1:
                        config[column[0]] = float(column[1])
                        parameters.config[column[0]] = float(column[1])
                    else:
                        config[column[0]] = int(column[1])
                        parameters.config[column[0]] = int(column[1])
    return config

def readline_binary(filename):
    line = numpy.fromfile(filename)
    n = len(line)/parameters.binfields
    line.shape = (n,parameters.binfields)
    return line

def plotline(freq,absorp,width,z):
    minfreq = min(freq)
    maxfreq = max(freq)
    xarray = numpy.arange(minfreq,maxfreq,parameters.max_resolve)
    yarray = numpy.zeros(len(xarray))
    for i in range(len(freq)):
        realwidth = width[i]/parameters.config['nu0']*freq[i]
        #print  absorp[i]
        yarray += absorp[i]/numpy.sqrt(2.*numpy.pi)/realwidth*numpy.exp(-0.5*((xarray-freq[i])/realwidth)**2)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xarray,yarray)
    plt.savefig('test.pdf')
    return yarray
    
    
