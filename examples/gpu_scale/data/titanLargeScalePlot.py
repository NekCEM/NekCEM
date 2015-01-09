import os
import matplotlib.pyplot as plt
import numpy as np
import sys

saveIn=sys.argv[1]
result=np.loadtxt('timingLargetitanGPUSort',dtype=[('No',np.int),('ele',np.int),('time',np.float)])
i=1

maxNo=16384
while i<=maxNo:
    plt.plot(result['No'][result['ele'] == i],result['time'][result['ele'] == i],'-o',label=str(i))
    i=8*i
plt.legend(loc='upper right')
i=16384
plt.plot(result['No'][result['ele'] == i],result['time'][result['ele'] == i],'-o',label=str(i))
plt.legend(loc='upper right')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Log(Number of GPUs)')
plt.ylabel('Log(Time (s))')
plt.title('Large Scale Titan Runs')
save='titanLargeScale'
if saveIn == '0':
    plt.savefig(save+'.pdf')
else:
    plt.show()
