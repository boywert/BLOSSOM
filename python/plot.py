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

def plotline(freq,absorp,width,output):
    minfreq = min(freq) - 20000.
    maxfreq = max(freq) + 20000.
    xarray = numpy.arange(minfreq,maxfreq,parameters.max_resolve)
    yarray = numpy.ones(len(xarray))

    for i in range(len(freq)):
        realwidth = width[i]/parameters.config['nu0']*freq[i]
        yarray *= 1.-absorp[i]*numpy.exp(-0.5*((xarray-freq[i])/realwidth)**2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xarray /= 1e6
    ax.plot(xarray,yarray)
    plt.savefig(output)
    return yarray

def plotline_sample(freq,absorp,width,output):
    minfreq = min(freq) - 20000.
    maxfreq = max(freq) + 20000.
    xarray = numpy.arange(minfreq,maxfreq,parameters.max_resolve)
    yarray = numpy.ones(len(xarray))

    for i in range(len(freq)):
        realwidth = width[i]/parameters.config['nu0']*freq[i]
        yarray *= 1.-absorp[i]*numpy.exp(-0.5*((xarray-freq[i])/realwidth)**2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(xarray)):
        xarray[i] -= freq[0]
    xarray /= 1e3
    ax.plot(xarray,yarray)
    plt.savefig(output)
    return yarray

def nu2Mpc(nu):
    H0 = 2.268545503707487E-018
    c = 29979000000.0
    nu0 = 1420405750.00000
    Mpc = 3.085677580000000E+024
    return c/H0*(1-(nu/nu0)**2.)/(1+(nu/nu0)**2.)/Mpc

def PrepFFT(freq,absorp,width):
    minfreq = min(freq)
    maxfreq = max(freq)
    nuarray = numpy.arange(minfreq,maxfreq,parameters.max_resolve)
    yarray = numpy.ones(len(nuarray))
    darray = numpy.ones(len(nuarray))
    print "N_freq",len(nuarray)
    for i in range(len(freq)):
        realwidth = width[i]/parameters.config['nu0']*freq[i]
        #yarray *= 1.-(1.-absorp[i])/numpy.sqrt(2.*numpy.pi)/realwidth*numpy.exp(-0.5*((nuarray-freq[i])/realwidth)**2)
        yarray *= 1.-(1.-absorp[i])*numpy.exp(-0.5*((nuarray-freq[i])/realwidth)**2)
    for i in range(len(yarray)):
        yarray[i] = max(0.,yarray[i])
        
    

    
