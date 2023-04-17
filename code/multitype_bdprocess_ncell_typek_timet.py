#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:39:24 2023

@author: mnichol3
"""



#Script to run a multitype birth-death process and number of cells of type k at given time

import numpy as np
import os
import random
# import matplotlib.pyplot as plt
from itertools import chain
import time
import pandas as pd
#import feather #using pyarrow instead
import pyarrow.feather as feather


#output filenames
#output_folder = "/exports/igmm/eddie/mnichol3-XDF/accumulate_nmutations/results/checking_variance/"
output_folder ="/home/mnichol3/ownCloud/accumulate_nmutations/simulation_python_output_local/fluctuation_assay_simouts/"



nruns = 10**2





seednum = 5
random.seed(seednum)

#we simulate until we hit nth type (but may as well save other data along way)
common_mutrate = 10**(-3)
ntypes = 3
fintime = 10 #add to string

shouldsave = True

output_filename = output_folder + 'popsizes' + '_seednum'+ \
                  str(seednum) + '_ntype' +str(ntypes) + '_nruns' + str(nruns) + \
                      '_fintime' +  str(fintime) + \
                       '_mrate' + str(common_mutrate) +'.feather' #(or .npy)

# seednum = 1
# random.seed(seednum)

#output_filename = output_folder + 'kstep4' + '_seednum'+ str(seednum) + '.npy'
#output_filename = output_folder + 'kstep4' + 'down_up' + '_seednum'+ str(seednum) + '.feather'
#output_filename = output_folder + 'twotype_popsizes' + '_seednum'+ str(seednum) + 'nrun'+ str(nruns) +'.feather'


#we simulate until we hit nth type (but may as well save other data along way)

#list_popsizes_overruns = [None] * nruns
array_popsizes_atfintime_overruns = np.empty([nruns,ntypes])
#set parameters for each type
#parameters when type 1 then down up
#birth_rates = np.array([1,1,1.4,1.03])
#death_rates = np.array([.3,1.5,.3,.1])

#parameters when type 1 then up dpwn 
birth_rates = np.array([1, 1, 1])
death_rates = np.array([0, 0, 0 ])

mut_rates = np.array([common_mutrate,common_mutrate,0])

total_rates_pertype = birth_rates + death_rates + mut_rates

#we wait until the first cell of type_to_hit exists. But if too long put in a max

# times_store_pop = np.array([.3,.6,.9,fintime])
#times_store_pop = np.array([.3,.6,.9,fintime])

#ntimes_store = len(times_store_pop)

actual_num_runs = 0
list_results_counter = 0
savecounter = 0

#nsavetimes = 10
#array_runtimestat_overruns = np.empty([nsavetimes,ntimes_store,6])


for run in range(nruns):

    print(run)
   
    pertype_cell_count = np.array([1,0,0]) #should be length of ntypes
    current_time = 0

    

    #passed_times = np.zeros(ntimes_store)
   # storing_times_passed = 0
    
    while ( (current_time<fintime) &  ( (pertype_cell_count[0]>0) | (pertype_cell_count[1]>0) | (pertype_cell_count[2]>0) )):
    
        #get next time event occurs
        eventrate = np.sum( total_rates_pertype*pertype_cell_count)
        deltat = np.random.exponential(scale = 1/eventrate)
        current_time = current_time + deltat
        
       
            
        if (current_time>= fintime):
            break
        
        #which type does new event occur in?
        type_event_rates = total_rates_pertype*pertype_cell_count
        
        multinomial_type_probs = np.array(type_event_rates )/eventrate
        logic_whichtype_event = np.random.multinomial(1,multinomial_type_probs)
        index_whichtype_event =  np.where( logic_whichtype_event  == 1)[0][0]
        
        
        #for chosen type get prob_birth, prob_death, prob_mut
        chosen_typebrate = birth_rates[index_whichtype_event]
        chosen_typedrate = death_rates[index_whichtype_event]
        chosen_typemut = mut_rates[index_whichtype_event]
        chosen_type_probdenom = chosen_typebrate + chosen_typedrate + chosen_typemut
                
        prob_birth =  chosen_typebrate/chosen_type_probdenom
        prob_death = chosen_typedrate/chosen_type_probdenom
        prob_mut=  chosen_typemut/chosen_type_probdenom
        
        multinomial_event_probs = [prob_birth,
                                           prob_death,
                                           prob_mut]
        logic_whichevent = np.random.multinomial(1, multinomial_event_probs)
        index_whichevent =  np.where( logic_whichevent   == 1)[0][0]
        
        
        #now update according to which event occured
        if index_whichevent == 0: #division
                    pertype_cell_count[index_whichtype_event] = pertype_cell_count[index_whichtype_event] + 1
                    
        elif index_whichevent == 1: #death
                    
                     pertype_cell_count[index_whichtype_event] = pertype_cell_count[index_whichtype_event] - 1 
                 
                    
        elif index_whichevent == 2: #mutation
                    
                     pertype_cell_count[index_whichtype_event+1]  = pertype_cell_count[index_whichtype_event+1] + 1
        
    array_popsizes_atfintime_overruns[run] =  pertype_cell_count


fintime_vec =  np.repeat([fintime],ntypes)
per_type_parameters = np.column_stack((birth_rates,death_rates,mut_rates,fintime_vec))

#first 3 columns are parameters, then sim results
params_sizes= np.concatenate([per_type_parameters,np.transpose(array_popsizes_atfintime_overruns)],axis=1)

params_sizes_pd  = pd.DataFrame(params_sizes)

if shouldsave == True:
    feather.write_feather(params_sizes_pd, output_filename)
