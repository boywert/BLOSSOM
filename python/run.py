import numpy
import matplotlib.pyplot as plt
import plot

configfile = '../inputs/config_20mpc_z8_ext1mpc'
plot.readconfig(configfile)

fig = plt.figure()


ax1 = fig.add_subplot(111)

linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RG/0.400/0/sout.1'

line = plot.readline_binary(linefile)
ak = numpy.histogram(line[:,4])
absorphist_y = ak[0] 
absorphist_x = ak[1]
absorphist_y *= 0
for j in range(4):
    for i in range(1,101):
        linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RG/0.400/%d/sout.%d'%(j,i)

        line = plot.readline_binary(linefile)
        a = numpy.histogram(line[:,4],absorphist_x)
        absorphist_y += a[0]
        print absorphist_y

#print len(absorphist_x),len(absorphist_y)
print sum(absorphist_y)
ax1.plot(absorphist_x[0:len(absorphist_y)],absorphist_y/float(sum(absorphist_y)),color = 'g')

ax1.set_xlim([absorphist_x[0],absorphist_x[len(absorphist_y)-1]])
ax1.set_yscale('log')
plt.savefig('hist_RG.pdf')



fig = plt.figure()
ax1 = fig.add_subplot(111)
linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RR/0.400/0/sout.1'

line = plot.readline_binary(linefile)
ak = numpy.histogram(line[:,4])
absorphist_y = ak[0] 
absorphist_x = ak[1]
absorphist_y *= 0
for j in range(4):
    for i in range(1,101):
        linefile = '/scratch/01937/cs390/outputs/cubepm_090610_14_5488_20Mpc/RESULT/ext1mpc/8.064/RR/0.400/%d/sout.%d'%(j,i)

        line = plot.readline_binary(linefile)
        a = numpy.histogram(line[:,4],absorphist_x)
        absorphist_y += a[0]
        print absorphist_y

print sum(absorphist_y)
ax1.plot(absorphist_x[0:len(absorphist_y)],absorphist_y/float(sum(absorphist_y)),color='b')



ax1.set_xlim([absorphist_x[0],absorphist_x[len(absorphist_y)-1]])
ax1.set_yscale('log')
plt.savefig('hist_RR.pdf')
#plot.plotline(line[:,3],line[:,4],line[:,5])
