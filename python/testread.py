#!/opt/apps/python/epd/7.2.2/bin/python
import numpy
#import matplotlib.pyplot as plt
import plot

configfile = '/home1/01937/cs390/BLOSSOM/inputs/config_20mpc_z8_ext1mpc'
plot.readconfig(configfile)





# for i in range(1,101):
#     print i
#     linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RG/0.000/0/sout.%d' % (i)
#     line = plot.readline_binary(linefile)
#     output = "python/000.%d.pdf" % (i)
#     for eachline in line:
#         print eachline
#     plot.plotline(line[:,3],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e5_point'
line = plot.readline_binary(linefile)
output = "python/sample_1e5_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e6_point'
line = plot.readline_binary(linefile)
output = "python/sample_1e6_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e7_point'
line = plot.readline_binary(linefile)
output = "python/sample_1e7_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e8_point'
line = plot.readline_binary(linefile)
output = "python/sample_1e8_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e5_extended'
line = plot.readline_binary(linefile)
output = "python/sample_1e5_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e6_extended'
line = plot.readline_binary(linefile)
output = "python/sample_1e6_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e7_extended'
line = plot.readline_binary(linefile)
output = "python/sample_1e7_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/Samples/Sample_1e8_extended'
line = plot.readline_binary(linefile)
output = "python/sample_1e8_p.pdf" 
for eachline in line:
    print eachline
    plot.plotline_sample(line[:,2],line[:,4],line[:,5],output)
