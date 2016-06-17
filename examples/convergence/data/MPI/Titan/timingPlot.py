import re
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

#Run with python timingPlot.py <save>
#if save is 0, it will save pngs.
#if save is 1, it will be interactive.

saveIn=sys.argv[1]
saveIn=int(saveIn)

plt.gcf().subplots_adjust(bottom=0.15)
#!!! Change font here
plt.rcParams.update({'font.size': 22})

#!!! Change file type here
fileExtension='.png'


for fn in os.listdir('.'):
    if os.path.isfile(fn):
        if(re.match('timing_comp.*P*_N7$',fn)):
            plt.clf()
            print fn
            result=np.loadtxt(fn,dtype=[('No',np.int),('ele',np.int),('time',np.float)])
            maxNo = max(result['No'])
            proc=fn.split("_")[3]
            comp=fn.split("_")[2]
            i=1
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*8*8*8,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('Time (s)')
            plt.xlabel('Number of Grid Points (N=7)')
            save='plots/'+fn+'lin'

            if saveIn == 0:
                plt.savefig(save+fileExtension)
            else:
                plt.show()

            i=1
            plt.clf()
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*8*8*8,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('Time (s)')

            plt.xscale('log')            
            plt.xlabel('Number of Grid Points (N=7)')

            save='plots/'+fn+'log'
            if saveIn == 0:
                plt.savefig(save+fileExtension)
            else:
                plt.show()

            plt.clf()
            i=1
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*8*8*8,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('Time (s)')
            
            plt.yscale('log')
            plt.xscale('log')          
            plt.xlabel('Number of Grid Points (N=7)')

            save='plots/'+fn+'loglog'
            if saveIn == 0:
                plt.savefig(save+fileExtension)
            else:
                plt.show()
