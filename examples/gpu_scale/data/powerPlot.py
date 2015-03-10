import re
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

#Run with python powerPlot.py <save>
#if save is 0, it will save pngs.
#if save is 1, it will be interactive.

saveIn=sys.argv[1]
saveIn=int(saveIn)

plt.gcf().subplots_adjust(bottom=0.15)
#!!! Change font here
plt.rcParams.update({'font.size': 20})

#!!! Change file type here
fileExtension='.png'

data = np.loadtxt('Cab1-noDate.csv')
data2 = np.loadtxt('Cab2-noDate.csv')
#calculate time and then start at 0
time = data[:,0]*3600 + data[:,1]*60 + data[:,2]
time = time - time[0]
#get GPU and CPU times
gpu = time[9380:9512]
gpu = gpu - gpu[0]
cpu = time[9645:9977]
cpu = cpu -cpu[0]
#Zero out the outer points to make a prettier plot
data[9645,3] = 0.0
data[9380,3] = 0.0
data2[9645,3] = 0.0
data2[9380,3] = 0.0
data[9976,3] = 0.0
data[9511,3] = 0.0
data2[9976,3] = 0.0
data2[9511,3] = 0.0
#plot data
plt.plot(cpu,data[9645:9977,3]+data2[9645:9977,3],label='CPU',color='blue',linewidth=2.5)
plt.plot(gpu,data[9380:9512,3]+data2[9380:9512,3],label='GPU',color='red',linewidth=2.5)
plt.title('NekCEM Power Consumption')
plt.xlabel('Time (s)')
plt.ylabel('Power (kW)')
plt.legend(loc='lower left')
plt.ylim(18,30)

save='plots/powerConsumption'
if saveIn == 0:
    plt.savefig(save+fileExtension)
else:
    plt.show()
