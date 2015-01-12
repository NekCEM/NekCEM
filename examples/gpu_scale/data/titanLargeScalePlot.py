import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from numpy.lib.recfunctions import append_fields


saveIn=sys.argv[1]
result=np.loadtxt('timing_comp_titan_large_GPU',dtype=[('No',np.int),('ele',np.int),('time',np.float)])
resultComm=np.loadtxt('timing_comm_titan_large_GPU',dtype=[('No',np.int),('ele',np.int),('time',np.float)])
#!!! Change font here
plt.rcParams.update({'font.size': 12})

#!!! Change file type here
fileExtension='.png'

i=1

maxNo=16384
while i<=maxNo:
    plt.plot(result['No'][result['ele'] == i],result['time'][result['ele'] == i],'-o',label=str(i))
    i=8*i

i=16384
plt.plot(result['No'][result['ele'] == i],result['time'][result['ele'] == i],'-o',label=str(i))
#plt.legend(loc='upper right')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of GPUs')
plt.ylabel('Time (s)')
plt.ylim([20,160])
plt.title('Large Scale Titan Runs')
save='plots/titanLargeScale'
if saveIn == '0':
    plt.savefig(save+fileExtension)
else:
    plt.show()


#--------------------------------------------
#
# Three horizontal lines weak scaling plot
#
#---------------------------------------------
plt.clf()
newResult = append_fields(result,'work',np.zeros(16))
newResultComm = append_fields(resultComm,'work',np.zeros(16))
i=1
while i<=maxNo:
    newResult['work'][newResult['ele'] == i] =  result['No'][result['ele'] == i]/float(i)
    newResultComm['work'][newResultComm['ele'] == i] =  resultComm['No'][resultComm['ele'] == i]/float(i)
    i=i*8

i=16384
newResult['work'][newResult['ele'] == i] =  result['No'][result['ele'] == i]/float(i)


i=0.25
while i<=1:
    plt.plot(newResult['No'][newResult['work'] == i],newResult['time'][newResult['work'] == i],'-o',label=str(i))
    i=i*2

#plt.legend(loc='upper right')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Number of GPUs')
plt.ylabel('Time (s)')
plt.title('Large Scale Titan Runs')
save='plots/titanLargeScale_horizontal'
if saveIn == '0':
    plt.savefig(save+fileExtension)
else:
    plt.show()

#--------------------------------------------
#
# Comp/comm weak scaling plot
#
#---------------------------------------------

plt.clf()

plt.plot(newResult['No'][newResult['work'] == 1],newResult['time'][newResult['work'] == 1],'-o',label='Computation')

plt.plot(newResultComm['No'][newResultComm['work'] == 1],newResultComm['time'][newResultComm['work'] == 1],'-o',label='Communication')

plt.legend(loc='lower right')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Number of GPUs')
plt.ylabel('Time (s)')
plt.title('Large Scale Titan Runs')
save='plots/titanLargeScale_comm'
if saveIn == '0':
    plt.savefig(save+fileExtension)
else:
    plt.show()
