#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 14:44:04 2020

@author: michaeldavidnicholson
"""

#Script to run a multitype birth-death process and get arrival times of type n cells
#script will need edited for birth, death, mutation parameters as required.
#output path needs edited also.
#saves output as 'feather' file, for ease loading to R for plotting.

import numpy as np
import random
import matplotlib.pyplot as plt
from itertools import chain
import time
import pandas as pd
from scipy import special
import feather


#output filenames

output_folder = "PATH_TO_PROJECTDIR" + "results/simout_arrival_times/"

#edit if want to save output
shouldsave = False


nsurvruns = 50


seednum = 2
random.seed(seednum)

#output_filename = output_folder + 'kstep4' + '_seednum'+ str(seednum) + '.npy'
output_filename = output_folder + 'kstep4' + 'down_up' + '_seednum'+ str(seednum) + '.feather'


#we simulate until we hit nth type (but may as well save other data along way)
type_to_hit = 4

result_hitting_times = np.zeros([nsurvruns,type_to_hit-1])

#set parameters for each type
#parameters when type 1 then down up
#birth_rates = np.array([1,1,1.4,1.03])
#death_rates = np.array([.3,1.5,.3,.1])

#parameters when type 1 then up dpwn 
birth_rates = np.array([1, 1.2 , 1.4, 1.6])
death_rates = np.array([0, 0 , 0,  0])
mut_rates = np.array([.1,.1,.1,.1])

total_rates_pertype = birth_rates + death_rates + mut_rates

#we wait until the first cell of type_to_hit exists. But if too long put in a max
max_time_sim_length = 30

actual_num_runs = 0

for run in range(nsurvruns):

    go_extinct = 0 
    
    if (run % 1 == 0):
        print(run)
        
    while (go_extinct == 0 ): 
            
        pertype_cell_count = np.array([1,0,0,0])
        current_time = 0
        

        #indicator if type i cell exists
        indicator_typeexists =  np.array([1,0,0,0])
        
        while (pertype_cell_count[type_to_hit-1])==0 & (pertype_cell_count[0]>0):
            
            #get next time event occurs
            eventrate = np.sum( total_rates_pertype*pertype_cell_count)
            deltat = np.random.exponential(scale = 1/eventrate)
            current_time = current_time + deltat
            
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
                         
                         
            
            #if first time populate new type, store that time
            if (pertype_cell_count[1] >0 ) & (indicator_typeexists[1] == 0):
                time_first_type2 = current_time
                indicator_typeexists[1] = 1 
                
            if (pertype_cell_count[2] >0 ) & (indicator_typeexists[2]  == 0):
                time_first_type3 = current_time
                indicator_typeexists[2] = 1 
                
            if (pertype_cell_count[3] >0 ) & ( indicator_typeexists[3] == 0):
                time_first_type4 = current_time
                indicator_typeexists[3] = 1 
        
        actual_num_runs  = actual_num_runs + 1
        if (pertype_cell_count[0] > 0):
               go_extinct = 1
    
    if type_to_hit == 2:
        result_hitting_times[run,0] = time_first_type2
    elif type_to_hit == 3:
        result_hitting_times[run,0] = time_first_type2
        result_hitting_times[run,1] = time_first_type3
    elif type_to_hit == 4:
        result_hitting_times[run,0] = time_first_type2
        result_hitting_times[run,1] = time_first_type3
        result_hitting_times[run,2] = time_first_type4
        
result_hitting_times_wzero = np.concatenate([np.zeros([nsurvruns,1]),result_hitting_times],axis = 1)
#np.save(output_filename,result_hitting_times)

result_hitting_times_wzero_t = np.transpose(result_hitting_times_wzero)
per_type_parameters = np.column_stack((birth_rates,death_rates,mut_rates))

#first 3 columns are parameters, then sim results
params_times = np.concatenate([per_type_parameters,result_hitting_times_wzero_t],axis=1)

params_times_pd  = pd.DataFrame(params_times)

if shouldsave == True:
    feather.write_dataframe(params_times_pd , output_filename)