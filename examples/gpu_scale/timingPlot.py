import re
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

save=sys.argv[0]

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
                plt.plot(result['time'][result['No'] == i],result['ele'][result['No'] == i], '-o',label=proc+str(i))
                i=2*i
            plt.legend(loc='lower right')
            plt.title('Timing Runs for '+comp+proc)
            plt.xlabel('Time (s)')
            plt.ylabel('No. of Elements')
            save=fn+'a'
            if save == 0:
                plt.savefig(save+'.png')
            else:
                plt.show()

