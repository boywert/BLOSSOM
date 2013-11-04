import numpy
#import matplotlib.pyplot as plt
import plot

configfile = '/home1/01937/cs390/BLOSSOM/inputs/config_20mpc_z8_ext1mpc'
plot.readconfig(configfile)


for i in range(1,101):
    print i
    linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RG/0.200/0/sout.%d' % (i)
    line = plot.readline_binary(linefile)
    for eachline in line:
        for lineel in eachline:
            if(lineel < 0.):
                print eachline

#plot.plotline(line[:,3],line[:,4],line[:,5],'0.000.pdf')
