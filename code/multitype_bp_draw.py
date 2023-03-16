# simulate a multytype process
# (i)->(i)(i) at rate ai, (i)->(i)(i+1) at mui, (i)->0 at rate bi
# start with a single (0);
# Produces Fig 1c in "Mutation accumulation in exponentially growing populations"

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

projectRoot = "/home/mnichol3/ownCloud/accumulate_nmutations/"

simOutDir = projectRoot + "results/simout/"
plotdir = projectRoot + "images/"

simOutFile = simOutDir+ "threeTypeOverTimePaperFig.csv"
plotFile = plotdir + "threeTypeOverTimePaperFig.pdf"

redoSims = 'F' #if True redo sims and save at simOutFile, else load sims from simOutFile.
saveplot = 'T' #if true save the putput figure


if redoSims == 'T':
    np.random.seed(10)
    
    numtypes = 3 # number of types
    a = np.zeros(numtypes, dtype=float)
    b = np.zeros(numtypes, dtype=float)
    mu = np.zeros(numtypes, dtype=float)
    n = np.empty(numtypes, dtype=int)
    a[0]= 1.1
    b[0]= 0.8
    mu[0]= 0.01
    a[1]= 1.0
    b[1]= 0.9
    mu[1]= 0.01
    a[2]= 1.1
    b[2]= 0.5
    mu[2]= 0.0
    
    stop = 35.0 # time when we stop the simulation
    
    # for plotting:
    maxstep = 100000 # max number of divisions we follow for each type
    steps = np.empty(numtypes, dtype=int) # step of number of cells for ech type
    times = np.empty((numtypes,maxstep+1), dtype=float) # times of jumps for each type
    sizes = np.empty((numtypes,maxstep+1), dtype=int) # number of cells after jump for each type
    
    if (mu[numtypes-1]>0): # just a check here
        print("last type should not mutate")
        exit()
    
    n.fill(0)
    while n[numtypes-1] == 0:
        sizes.fill(0)
        steps.fill(0)
        n.fill(0)
        times.fill(np.nan) # I need these to avoid funny plots at the end
        for i in np.arange(0,numtypes):
            times[i,0] = 0
        n[0] = 1 # initially we have one wild-type cell
        sizes[0,0] = 1 # first type starts at 1
    
        totalrate = 0
        for i in np.arange(0,numtypes):
            totalrate += n[i]*(a[i]+b[i]+mu[i])
    
        time = np.random.exponential(1./totalrate)
        while time<stop and totalrate>0:
            ra = np.random.random() * totalrate # store the random number
            nowtype=-1 # pick the type which will change
            while(ra>0):
                nowtype += 1
                ra -= n[nowtype]*(a[nowtype]+b[nowtype]+mu[nowtype])
            # nowtype will act
            ra += n[nowtype]*a[nowtype]
            if (ra>0): # i divides
                n[nowtype] += 1
            else:
                ra += n[nowtype]*b[nowtype]
                if (ra>0): # nowtype dies
                    n[nowtype] -= 1
                else: # nowtype mutatates
                    nowtype += 1
                    n[nowtype] += 1
    
            # For plotting:
            if steps[nowtype]<maxstep:
                steps[nowtype] += 1
                times[nowtype, steps[nowtype]] = time
                sizes[nowtype, steps[nowtype]] = n[nowtype]
    
            totalrate = 0
            for i in np.arange(0,numtypes):
                totalrate += n[i]*(a[i]+b[i]+mu[i])
            if totalrate>0:
                time += np.random.exponential(1./totalrate)
    
    
    
    #now fit W_i* exp(lambda_i*t)
    def get_random_amp(reltimes,relsizes,relgrate):
        # reltimes = times[1]
        # relsizes = sizes[1]
        # relgrate = a[0]-b[0]
        #strip nan
        reltimes = [reltimes[x] for x in range(len(reltimes)) if str(reltimes[x])!='nan']
        relsizes = relsizes[0:(len(reltimes)-1)]
        fintime = reltimes[-1]
        indexhit =next(x[0] for x in enumerate(relsizes) if x[1] >= 1*(10**2))
        
        #amplitude given by (y2-y1)/(exp(l*t2)-exp(l*t1))
        W = (relsizes[-1]-relsizes[indexhit])/(np.exp(relgrate*reltimes[-1]) - np.exp(relgrate*reltimes[indexhit]))
        
        return(W)
    
    grates = a-b
    running_max = [grates[0],grates[0],grates[2]]
    
    amplitudes = [get_random_amp(times[i], sizes[i], running_max[i]) for i in range(len(times))]
    
    asymp_sizes = [amplitudes[i]*np.exp(running_max[i]*times[i]) for i in range(len(times))]
    
    
    times_sizes= np.concatenate([times,sizes,asymp_sizes],axis=0)
    
    times_sizes_pd0  = pd.DataFrame(times_sizes)
    times_sizes_pd1 = times_sizes_pd0.transpose()
    times_sizes_pd1.columns =['times1','times2','times3','sizes1','sizes2','sizes3','asymsizes1', 'asymsizes2','asymsizes3']
    # times_sizes_pd1.to_csv(simOutFile)
    
else:
        
    times_sizes_pd1 = pd.read_csv(simOutFile)

if saveplot == 'T':
    
    Path(plotdir).mkdir(parents=True, exist_ok=True) #create image directory if doesn't exist yet
    #print(steps, times, sizes)
    numtypes = int( (len(times_sizes_pd1.columns)-1)/3)
    colors = ['b','r','y']
    plt.xlim([0, 35])
    plt.ylim([0.8, 10000])
    plt.yscale("log")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$Z_n(t)$')
    
    for i in np.arange(1,numtypes+1):
        plt.step( times_sizes_pd1.iloc[:,[i]], times_sizes_pd1.iloc[:,[3+i]], where='post',color = colors[i-1])
        plt.plot(times_sizes_pd1.iloc[:,[i]], times_sizes_pd1.iloc[:,[6+i]], color = colors[i-1],linestyle = 'dashed')
     
    
    
    figure = plt.gcf()
    figure.set_size_inches(8.5, 3)
    plt.savefig(plotFile )
    plt.show()