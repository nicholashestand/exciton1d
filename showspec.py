#Script to plot spectrum using matplotlib
import sys
import matplotlib.pyplot as plt
import numpy as np

task_title = sys.argv[1]
specf = task_title + '_ab.csv'
dispf = task_title + '_disp.csv'
spec = np.loadtxt(specf,delimiter=',',skiprows=1)
disp = np.loadtxt(dispf,delimiter=',',skiprows=1)

#normalize dispersion data to pi
disp[:,0] = 2*disp[:,0]/len(disp)

#plot the spectrum with matplotlib
plt.subplot(2,1,1)
plt.plot( spec[:,0], spec[:,1], c='r' )
plt.xlabel(r'Energy (h$\omega$)')
plt.ylabel('Absorption (a.u.)')
plt.title('Absorption Spectrum')

#plot the dispersion with matplotlib
plt.subplot(2,1,2)
plt.plot( disp[:,0], disp[:,1], c='b' )
plt.xlabel(r'k ($\pi$)')
plt.ylabel(r'Energy (h$\omega$)')
plt.title('Lowest Exciton Dispersion')
plt.tight_layout()
plt.show()
plt.close()
