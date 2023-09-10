# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:35:23 2020

@author: maryam
"""
#Homework 4 - Evolutionary Strategies 
"""
To do this you need to: 
    1- encode the problem as two real variables
    2- select values for μ and lambda
    3- select a standard deviation (the same for both variables) and change it 
    using the 1/5 success rule
    4- randomly select an initial population
    5- select a stopping criterion. 
    --> Do not perform recombination.
"""

#%% #import libraries 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random 
np.random.seed(256)

#%% # initialize mu, lamda, and the ratio 
ratio = 8
#number of parents
mu = 10 
#number of children
lamda = mu * ratio

sigma = 1.0
alpha = 0.84

#stopping criteria
runs = 100

#%%
x = range (-6, 6)
y = range (-6, 6)

lower_b = -6
upper_b = 6 
#%% #Himmelblau function”
def himmelblau(x,y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

#%%
##################################################################
    #without recombination 
###################################################################
    
solution = {}

selected_parent = 0
child = {}
selected_parent_dic = {}

all_parent_child = []
df_data = pd.DataFrame(columns = ['g', 'j','r', 'fitness_child','child[0]', 'child[2]'])


#%%

row_counter = 0 
for g in range(runs):
      
    #def generate_parent():
    # Set a length of the list to k
    parent = np.zeros((mu, 4), dtype=np.float64)
    

    
    for j in range (mu):
        child[j] = {}
        selected_parent_dic[j] = {}
        for i in range(mu):
        
            rand_x = np.random.uniform(lower_b,upper_b)
            rand_y = np.random.uniform(lower_b,upper_b)
            
            parent[i][0] = rand_x
            parent[i][1] = sigma
            parent[i][2] = rand_y
            parent[i][3] = sigma
    #        print(parent[i])
    
        selected_parent = parent[random.randint(0, 9)]
        selected_parent_dic[j] = selected_parent
    
    #    list_p = []
    #    list_c = []
        #child = [None]*4
         
        for r in range (ratio):
            child[j][r]={}
            child[j][r][0] = selected_parent[0] + np.random.normal(0, sigma)
            child[j][r][1] = sigma
            child[j][r][2] = selected_parent[2] + np.random.normal(0, sigma)
            child[j][r][3] = sigma
#            print('=================================')
#            print(child[j][r])
#            print('++++++++++++++++++++++++++++++++++++')
#            print(child[j][r][0], child[j][r][2])
#            print('**********************************')
    
            counter = 0
            fitness_parent = himmelblau(selected_parent[0],selected_parent[2])
            fitness_child = himmelblau (child[j][r][0], child[j][r][2])
            
            

            df_data.at[row_counter, 'fitness_child'] = fitness_child
            df_data.at[row_counter, 'r'] = r
            df_data.at[row_counter, 'j'] = j
            df_data.at[row_counter, 'child[0]'] = child[j][r][0]
            df_data.at[row_counter, 'child[2]'] = child[j][r][2]
            df_data.at[row_counter, 'g'] = g


            all_parent_child.append(fitness_child)
            
            row_counter += 1
            
            if fitness_child > fitness_parent:
                counter += 1
                
                success = counter / lamda
                if success > 1/5:
                    sigma = (1/alpha)*sigma
                elif success< 1/5:
                    sigma = alpha *sigma   
                            
#%%
df_data['fitness_child'].min()

#%%
a = df_data.sort_values(by=['fitness_child'], ascending=True).iloc[0:10]
a

#%%

################################################################################
# Add Global Intermediate
################################################################################

parents = np.zeros((5, 4), dtype=np.float64)
sigmax = 1.0 
sigmay = 2.0
#%%
def globl_intermediate():
#
    for p in range (5):

        rand_x = np.random.uniform(lower_b,upper_b)
        rand_y = np.random.uniform(lower_b,upper_b)
        
        parents[p][0] = rand_x
        parents[p][1] = sigmax
        parents[p][2] = rand_y
        parents[p][3] = sigmay
         
#
# Add global intermediate recombination 
#    first_parent = parent[0]
#    second_parent = [parent[1][0], parent[2][1], parent[3][2], parent[4][3]]
    
    new_x = (parents[0][0] + parents[1][0])/2
    new_sigx = (parents[0][1] + parents[2][1])/2
    new_y = (parents[0][2] + parents[3][2])/2
    new_sigy = (parents[0][3] + parents[4][3])/2
    
    recombined_parent = [new_x, new_sigx, new_y, new_sigy]
 #  
    return recombined_parent
    
    
#%%    

solution = {}

selected_parent = 0
child = {}
selected_parent_dic = {}

all_parent_child = []

df_data = pd.DataFrame(columns = ['g', 'j','r', 'fitness_child','child[0]', 'child[2]'])


#%%

row_counter = 0
for g in range(runs):
#def generate_parent():
# Set a length of the list to k
    parent = np.zeros((mu, 4), dtype=np.float64)
    
    
    for j in range (mu):
        child[j] = {}
        selected_parent_dic[j] = {}
    
        selected_parent = globl_intermediate()
        selected_parent_dic[j] = selected_parent
    
    #    list_p = []
    #    list_c = []
        #child = [None]*4
         
        for r in range (ratio):
            child[j][r]={}
            child[j][r][0] = selected_parent[0] + np.random.normal(0, sigma)
            child[j][r][1] = sigmax
            child[j][r][2] = selected_parent[2] + np.random.normal(0, sigma)
            child[j][r][3] = sigmay
            print('=================================')
            print(child[j][r])
            print('++++++++++++++++++++++++++++++++++++')
            print(child[j][r][0], child[j][r][2])
            print('**********************************')
    
            counter = 0
            fitness_parent = himmelblau(selected_parent[0],selected_parent[2])
            fitness_child = himmelblau (child[j][r][0], child[j][r][2])
            
            df_data.at[row_counter, 'fitness_child'] = fitness_child
            df_data.at[row_counter, 'r'] = r
            df_data.at[row_counter, 'j'] = j
            df_data.at[row_counter, 'child[0]'] = child[j][r][0]
            df_data.at[row_counter, 'child[2]'] = child[j][r][2]
            df_data.at[row_counter, 'g'] = g
            
            all_parent_child.append(fitness_child)
            
            row_counter += 1
            
            if fitness_child > fitness_parent:
                counter += 1
                
                success = counter / lamda
                if success > 1/5:
                    sigma = (1/alpha)*sigma
                elif success< 1/5:
                    sigma = alpha *sigma   
    
#%%     
df_data['fitness_child'].min()

#%%
a = df_data.sort_values(by=['fitness_child'], ascending=True).iloc[0:10]
a

#%%

################################################################################
# Add Dual Discrete
################################################################################
parents = np.zeros((2, 4), dtype=np.float64)
new_parent = np.zeros((1, 4), dtype=np.float64)
sigmax = 1.0 
sigmay = 2.0

#%%
def discrete_dual():
    
    for d in range (2):

        rand_x = np.random.uniform(lower_b,upper_b)
        rand_y = np.random.uniform(lower_b,upper_b)
        
        parents[d][0] = rand_x
        parents[d][1] = sigmax
        parents[d][2] = rand_y
        parents[d][3] = sigmay

    k1 = random.random()
    if k1 < 0.5:
        x = parents[0][0]
    else: 
        x = parents[1][0]
    
    k2 = random.random()
    if k2 < 0.5:
        sig_x = parents[0][1]
    else: 
        sig_x = parents[1][1]
        
    k3 = random.random()
    if k3 < 0.5:
        y = parents[0][2]
    else: 
        y = parents[1][2]
        
    k4 = random.random()
    if k4 < 0.5:
        sig_y = parents[0][3]
    else: 
        sig_y = parents[1][3]
       
    new_parent = [x, sig_x, y, sig_y]
    
    return  new_parent
        
#%%
solution = {}

selected_parent = 0
child = {}
selected_parent_dic = {}

all_parent_child = []

df_data = pd.DataFrame(columns = ['g', 'j','r', 'fitness_child','child[0]', 'child[2]'])


#%%

row_counter = 0
for g in range(runs):
#def generate_parent():
# Set a length of the list to k
    parent = np.zeros((mu, 4), dtype=np.float64)
    
    
    for j in range (mu):
        child[j] = {}
        selected_parent_dic[j] = {}
    
        selected_parent = discrete_dual()
        selected_parent_dic[j] = selected_parent
    
    #    list_p = []
    #    list_c = []
        #child = [None]*4
         
        for r in range (ratio):
            child[j][r]={}
            child[j][r][0] = selected_parent[0] + np.random.normal(0, sigma)
            child[j][r][1] = sigmax
            child[j][r][2] = selected_parent[2] + np.random.normal(0, sigma)
            child[j][r][3] = sigmay
            print('=================================')
            print(child[j][r])
            print('++++++++++++++++++++++++++++++++++++')
            print(child[j][r][0], child[j][r][2])
            print('**********************************')
    
            counter = 0
            fitness_parent = himmelblau(selected_parent[0],selected_parent[2])
            fitness_child = himmelblau (child[j][r][0], child[j][r][2])
            
            df_data.at[row_counter, 'fitness_child'] = fitness_child
            df_data.at[row_counter, 'r'] = r
            df_data.at[row_counter, 'j'] = j
            df_data.at[row_counter, 'child[0]'] = child[j][r][0]
            df_data.at[row_counter, 'child[2]'] = child[j][r][2]
            df_data.at[row_counter, 'g'] = g
            
            all_parent_child.append(fitness_child)
            
            row_counter += 1
            
            if fitness_child > fitness_parent:
                counter += 1
                
                success = counter / lamda
                if success > 1/5:
                    sigma = (1/alpha)*sigma
                elif success< 1/5:
                    sigma = alpha *sigma   
    
#%%     
df_data['fitness_child'].min()

#%%
a = df_data.sort_values(by=['fitness_child'], ascending=True).iloc[0:10]
a










