import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial2
import os

pathToMLdata = '/home/mnichol3/ownCloud/accumulate_nmutations/results/MLfigure_data'
os.chdir(pathToMLdata)

cadv = "#000000"
cless = "#0072B2"
cdis = "#009E73"

plt.xlabel(r'$V_n/\omega_n$')
plt.ylabel(r'density')
plt.xscale("log")
plt.yscale("log")
plt.xlim(0.001,200)
plt.ylim(0.0001,10)
#sim6a[0,1] = np.nan
#plt.plot(sim6a[:,0]/scale2a,sim6a[:,1]*scale2a,'-', label='simulation')

#2 species, 7.000000 maxi, and 100000 runs
#first species: d-b-m: 0.000000 1.000000 0.010000
#second species: d-b-m: 0.000000 2.000000 0.000000
# scale = 296.731
# with open('sim_2typeb.txt') as f:
# da = np.array(x2)
# digitized = np.digitize(da[:,0], bins)
# bin_means = np.array([da[digitized == i,1].mean() for i in range(1, len(bins))])
# #plt.plot(da[:,0],da[:,1],'x', label='simulation')
# plt.plot(bin_middles/scale, bin_means*scale,'o', label='simulation')


# 2 species, 5.000000 maxi, and 10000 runs
# first species: d-b-m: 0.200000 1.200000 0.010000
# second species: d-b-m: 0.300000 4.300000 0.000000
bins = [0,1,2,3,4,5,7,10,20,40,80,160,300,1000,3000,10000,30000,100000]
bin_middles = np.array([np.sqrt(bins[i]*(bins[i+1]-1)) for i in range(0, len(bins)-1)])
scale = 12.3256
with open('sim_2typee.txt') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
da = np.array(x2)
digitized = np.digitize(da[:,0], bins)
bin_means = np.array([da[digitized == i,1].mean() for i in range(1, len(bins))])
plt.plot(bin_middles/scale, bin_means*scale,'o', color=cadv, label=r'sim 1, $\lambda_1/\delta_2 = 0.25$')


#2 species, 7.000000 maxi, and 100000 runs
#first species: d-b-m: 0.200000 1.200000 0.010000
#second species: d-b-m: 0.300000 2.300000 0.000000
bins = [0,1,2,3,4,5,7,10,20,40,80,160,300,1000,3000,10000,30000,100000]
bin_middles = np.array([np.sqrt(bins[i]*(bins[i+1]-1)) for i in range(0, len(bins)-1)])
scale = 371.558
with open('sim_2typec.txt') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
da = np.array(x2)
digitized = np.digitize(da[:,0], bins)
bin_means = np.array([da[digitized == i,1].mean() for i in range(1, len(bins))])
plt.plot(bin_middles/scale, bin_means*scale,'o', color=cless, label=r'sim 2, $\lambda_1/\delta_2 = 0.5$')


#2 species, 12.000000 maxi, and 10000 runs
#first species: d-b-m: 0.200000 1.200000 0.010000
#second species: d-b-m: 0.300000 1.000000 0.000000
bins = [0,1,2,3,4,5,6,7,10,20,40,80,160,300,1000,2000,10000,15000,30000,60000,80000, 100000]
bin_middles = np.array([np.sqrt(bins[i]*(bins[i+1]-1)) for i in range(0, len(bins)-1)])
scale = 6510.19
with open('sim_2typed.txt') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
da = np.array(x2)
digitized = np.digitize(da[:,0], bins)
bin_means = np.array([da[digitized == i,1].mean() for i in range(1, len(bins))])
plt.plot(bin_middles/scale, bin_means*scale,'o', color=cdis, label=r'sim 3, $\lambda_1/\delta_2 = 1.0$')




# 3 species, 12.000000 maxi, and 10000 runs
# first species: d-b-m: 0.200000 1.200000 0.010000
# second species: d-b-m: 0.300000 1.000000 0.001000
bins = [0,1,2,3,6,10,20,40,80,160,300,1000,3000,10000,30000,100000]
bin_middles = np.array([np.sqrt(bins[i]*(bins[i+1]-1)) for i in range(0, len(bins)-1)])
scale = 87.1457
with open('sim_2typeg.txt') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
da = np.array(x2)
digitized = np.digitize(da[:,0], bins)
bin_means = np.array([da[digitized == i,1].mean() for i in range(1, len(bins))])
plt.plot(bin_middles/scale, bin_means*scale,'s', color=cless, label=r'sim 4, $\lambda_1/\delta_3 = 0.5$')


# --------------------------------------
# Theory curves:
with open('ML0.25.dat') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
theo5 = np.array(x2)
plt.plot(theo5[:,0],theo5[:,1],'-', color=cadv)

with open('ML0.5.dat') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
theo5 = np.array(x2)
plt.plot(theo5[:,0],theo5[:,1],'-', color=cless)

with open('ML1.dat') as f:
   x2 = [[float(x) for x in line.split()] for line in f]
theo5 = np.array(x2)
plt.plot(theo5[:,0],theo5[:,1],'-', color=cdis)


plt.legend(loc='lower left')
#plt.savefig('fig_Mittag.pdf')
plt.show()
