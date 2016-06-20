import matplotlib
matplotlib.use('TkAgg')
import sys
from matplotlib import gridspec
from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 16})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lw=2
linestyle = ["g--","k:","r-","c-.","g:","k-","r-.","c--","g-","k-.","r--","c-"]

#Run with python titan_plot.py

#import data
data = np.loadtxt("timing_comp_Titan_MPI")


fig, ax = plt.subplots()

plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)

# add some text for labels
ax.set_ylabel('Time (s)',fontsize=22)
ax.set_xlabel('Data Size',fontsize=22)
# title
ax.set_title('Time vs Data Size',fontsize=22)


j=0
for i in [2,4,8,10,12,14,16,18]:
    xdata = data[data[:,1]==i,0]**3*data[data[:,1]==i,1]**3
    ydata = data[data[:,1]==i,2]
    plt.plot(xdata,ydata,linestyle[j],label='N='+str(i),linewidth=lw)
    j=j+1

# We first sort the data to allow us to plot a line in the correct order
# Ordered lines
# xdata = data[:,0]**3*data[:,1]**3
# indices = np.argsort(xdata)

# plt.plot(xdata[indices],data[indices,2],label='Data',linewidth=lw)

#Set logscale
plt.yscale('log')
plt.xscale('log')
#Show legend
plt.legend(loc="upper left")
save = 'titan_timing'

plt.savefig(save+'.pdf')

#Now plot error
plt.clf()
fig, ax = plt.subplots()

plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)

# add some text for labels
ax.set_ylabel('Error',fontsize=22)
ax.set_xlabel('Data Size',fontsize=22)
# title
ax.set_title('Error vs Data Size',fontsize=22)

j=0
for i in [2,4,8,10,12,14,16,18]:
    xdata = data[data[:,1]==i,0]**3*data[data[:,1]==i,1]**3
    ydata = data[data[:,1]==i,3]
    plt.plot(xdata,ydata,linestyle[j],label='N='+str(i),linewidth=lw)
    j=j+1

plt.yscale('log')
plt.xscale('log')
plt.legend(loc="lower left")

save = 'titan_error'

plt.savefig(save+'.pdf')

