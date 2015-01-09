import re
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

#Run with python timingPlot.py <save>
#if save is 0, it will save pngs.
#if save is 1, it will be interactive.

saveIn=sys.argv[1]
plt.rcParams.update({'font.size': 22})
plt.gcf().subplots_adjust(bottom=0.15)

for fn in os.listdir('.'):
    if os.path.isfile(fn):
        if(re.match('timing.*Sort$',fn)):
            plt.clf()
            result=np.loadtxt(fn,dtype=[('No',np.int),('ele',np.int),('time',np.float)])
            maxNo = max(result['No'])
            proc=fn[-7:-4]
            comp=fn[6:-7]
            i=1
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*15*15*15,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('Time (s)')
            plt.xlabel('Number of Grid Points')
            save=fn+'lin'
            if saveIn == '0':
                plt.savefig(save+'.pdf')
            else:
                plt.show()

            i=1
            plt.clf()
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*15*15*15,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('Time (s)')

            plt.xscale('log')            
            plt.xlabel('log(Number of Grid Points)')

            save=fn+'log'
            if saveIn == '0':
                plt.savefig(save+'.pdf')
            else:
                plt.show()

            plt.clf()
            i=1
            while i<=maxNo:
                plt.plot(result['ele'][result['No'] == i]*15*15*15,result['time'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i

            plt.legend(loc='upper left',prop={'size':20})
            plt.title('Timing Runs for '+comp+proc)
            plt.ylabel('log(Time (s))')
            
            plt.yscale('log')
            plt.xscale('log')          
            plt.xlabel('log(Number of Grid Points)')

            save=fn+'loglog'
            if saveIn == '0':
                plt.savefig(save+'.pdf')
            else:
                plt.show()
