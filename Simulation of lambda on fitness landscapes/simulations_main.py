#!/usr/bin/env python
# coding: utf-8

# In[2]:


"""To simulate lambda evolution on ancestor host (606) and resistant host (Ecc4) fitness landscapes"""

#To get performance bottlenecks and python code as heatmap. Make sure to install py-heat first ("pip install py-heat")
# load_ext heat 
# %%heat

import igraph

import numpy as np
from numpy.random import binomial as nbinom
from numpy.random import multinomial as nmultinom
from numpy.random import choice as nchoice
from numpy.random import gamma as rgamma
from scipy.stats import gamma as pgamma
from random import choices
import random
import math

from timeit import default_timer as timer

from numpy import nonzero as nnonzero
from six.moves import range as srange
import csv
import pandas as pd
from matplotlib import pyplot as plt

import plotly.plotly as py
import plotly.figure_factory as ff

from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('pdf', 'svg')



def make_genotype_space(num_loci):
  genotype_space = igraph.Graph(directed=False, 
                      graph_attrs={'num_loci':num_loci, 'generation':0})

  genotype_space.add_vertices(2**num_loci)  # creates a graph with 2**num_loci vertices

  edges_to_add = []  # list which defines nearest neighbours
  for i in srange(0, 2**num_loci):  
    for j in srange(i, 2**num_loci):
      if bin(i^j).count('1') == 1: # ^ is bitwise XOR comparison to find nearest neighbors. NOTE the genotype ID to binary string for mutation is given by int_to_bin function which is mirror image of this n=binary representation. The neighbour connections does not change so we can use this here.
        edges_to_add.append( (i,j) )

  genotype_space.add_edges(edges_to_add)
  distance_to_optimum = [num_loci - bin(i^(2**num_loci - 1)).count('1') for i in srange(2**num_loci)] # no. of mutations in the genotype
  genotype_space.vs['distance'] = distance_to_optimum
  genotype_space.vs['label'] = [i for i in srange(2**num_loci)]  # srange goes from 0 to 2**num_loci - 1
  genotype_space.vs['fitness'] = [None for i in srange(2**num_loci)]
  genotype_space.vs['abundance'] = [0 for i in srange(2**num_loci)]
  genotype_space.vs['death_rate'] = [0 for i in srange(2**num_loci)]  

  return genotype_space

def nearest_neighbour_imputation(node, network):
    neighbour_fitnesses = np.array(network.vs.select(network.neighbors(node))['fitness'],dtype=float)
    non_nan_neighbour_fitnesses = neighbour_fitnesses[~np.isnan(neighbour_fitnesses)]  # removes all the nan fitness values
    if non_nan_neighbour_fitnesses.size > 0:
        return np.mean(non_nan_neighbour_fitnesses,dtype=float)
    else:
        return None
 

def add_residual_noise(switch, residual_matrix):
    if switch:
        residual_matrix = residual_matrix[~np.isnan(residual_matrix)]
        residual_noise = np.asscalar(random.choice(residual_matrix))
        return residual_noise
    else:
        return 0
    
    
    
def randomly(seq):
    shuffled = list(seq)
    random.shuffle(shuffled)
    return iter(shuffled)

def load_landscape(landscape_dataframe, wildtype_fitness, fitness_correction, imputation_function): # removinng residual matrix-> , residual_matrix): # landscape_dataframe is input of selection rates from a csv file with 2 columns- Genotype (0-1023) and SelectionRate
  max_genotype = np.max(landscape_dataframe.Genotype) # getting the maximum genotype ID or total no. of genotypes present
  genotype_network = make_genotype_space(int(np.log2(max_genotype+1))) # Creating igraph (genotype network) with input as num_loci (or no. of mutation loci)
  for id, row in landscape_dataframe.iterrows(): 
    
    selectionRates  = np.array([row.SelectionRate1, row.SelectionRate2, row.SelectionRate3, row.SelectionRate4])
    selectionRates = selectionRates[~np.isnan(selectionRates)]  # removing all the nan values from the selectionRates
    avg_selection_rate = np.nanmean(choices(selectionRates,k=len(selectionRates)),dtype=float)  # getting the average selection rate
    genotype_network.vs[int(row.Genotype)]['fitness'] = 0.5*0.25 * (wildtype_fitness + avg_selection_rate + fitness_correction )  # This is fitness per hr as the the selection rates are for 4 hrs

  no_fitness_nodes = genotype_network.vs.select(fitness_eq=None)

  while np.array(no_fitness_nodes).size > 0:
        no_fitness_nodes = genotype_network.vs.select(fitness_eq=None)
        for node in randomly(no_fitness_nodes):
            node['fitness'] = imputation_function(node, genotype_network)
#             if node['fitness']is not None:
#                 node['fitness'] = node['fitness'] + add_residual_noise(False, residual_matrix)
#     print(type(node))
#     print(node['fitness'])
#   for node in genotype_network.vs:
#     node['death_rate'] = death_rate_input.DeathRate[node['distance']]
#     genotype_space.vs['death_rate'] = [death_rate_input[node[]]] 
  return genotype_network    

# This will save the genotype and the selection rates (after imputing for the missing genotypes) in a csv file.
# Note that it's modified so that only selection rates from each run are saved and not the 1st column of genotype. A existing csv file with Genotype column is needed
def save_landscape(landscape, run_number, filename):
  dict_array = []
  wt_node = landscape.vs[0];
  for node in landscape.vs:
    dict_array.append({'Genotype':node['label'], 'SelectionRate':node['fitness'] - wt_node['fitness']})  # need to put fitness of WT instead of 1
  df = pd.DataFrame(dict_array)
  if run_number == 0:
    df.to_csv(filename,  index=False)
  else:  
    df_old = pd.read_csv(filename)
    df_old['SelectionRate'+ str(run_number+1)] = df['SelectionRate'];
    df_old.to_csv(filename,  index=False)

# This will pick a neighbour randomly to assign it to the mutant genotype
def get_mutant(population, parent_genotype):
  neighbors = population.neighbors(parent_genotype)
  mutant_id = nchoice(neighbors)
  return mutant_id

#this is Wright-Fisher dynamics implemented with binomial sampling, 
#which is a multinomial when you consider all the individual sub populations
def reproduce(landscape):
    population_size = landscape['population_size']
    fitnesses = np.array(landscape.vs['fitness'])
    abundances = np.array(landscape.vs['abundance'])
    weighted_fit = np.exp(fitnesses) * abundances

#     ab_weighted_fitness = fitnesses * (abundances/float(population_size))


    landscape.vs['abundance'] = nmultinom(n=population_size,
                                        pvals=weighted_fit / weighted_fit.sum(),
                                        size=1)[0]
#     landscape['mean_fitness'] = ab_weighted_fitness.sum()
    landscape['generation'] += 1
    return landscape

#At the given mutation rate, assign abundance to new genetic neighbors
def mutate(population, mutation_rate):
  """Mutate individuals in the population"""
  assert mutation_rate >= 0 and mutation_rate <= 1
#Use number_replication_genotypes if using death rate and above 3 lines
#   number_replications_genoptypes =  np.array( population.vs['abundance'] * (np.exp(np.array(population.vs['death_rate']) + np.array(population.vs['fitness'])) - 1) )
#   number_replications_genoptypes = number_replications_genoptypes.clip(min=0)
#   number_replications_genoptypes = number_replications_genoptypes.astype(int)

#num-mutants follow binomial distribution for each genotype ID of existing mutants; length of num_mutants is 1024
  num_mutants = nbinom(n=population.vs['abundance'], 
                       p=1 - np.exp(-mutation_rate))     
  population.vs['abundance'] = population.vs['abundance'] - num_mutants  
  for genotype_id in nnonzero(num_mutants)[0]:
    for i in srange(num_mutants[genotype_id]):
      new_mutant = get_mutant(population, genotype_id)
      population.vs[genotype_id]['abundance'] -= 1
      population.vs[new_mutant]['abundance'] += 1  
  return(num_mutants)

def plot_population(population, filename=None):  # Why filename=None here. Does it get overwritten when you provide the function with a filename
  layout = [[v['distance'], -v['fitness']] for v in population.vs]
  population.vs['size'] = [np.log(v['abundance']+0.01)*5 for v in population.vs]  # to show the relative abundance of genotypes by size of red circles 

  igraph.plot(population, target=filename,layout=layout)

#This is where the evolution actually happens...
#Really simple, just reproduce and mutate stuff
def do_timestep(population, mut_rate, num_generations, OmpF_ID): #, plot_population, t): # the input for num_generations will overwrite num_generations=1 
  OmpF_abundance = np.zeros([num_generations, len(OmpF_ID)]) 
  for gen in srange(num_generations):
#     plot_population(population, 'phage_{0}.png'.format(t+gen))    
    reproduce(population)
    num_mutants_do_timestep = mutate(population, mut_rate)  # num_mutants_do_timestep is sum of total mutations using multinomial function
    OmpF_abundance[gen,:] = [population.vs[x]['abundance'] for x in OmpF_ID]
  return(OmpF_abundance) # num_mutants_do_timestep

def int_to_bin(int_phenotype, num_loci):
  return bin(int_phenotype)[2:].zfill(num_loci)[::-1]

def bin_to_int(str_binary_phenotype):
  total = 0
  for i in srange(len(str_binary_phenotype)):
    if str_binary_phenotype[i] == '1':
      total += 2**i
  return total

def check_ompF(int_id, num_loci):
    z=0
    bin_id = int_to_bin(int_id, num_loci);
    ompF_required_bin = '0000101010'

    ompF_array = np.array([int(i) for i in ompF_required_bin])
    bin_array = np.array([int(i) for i in bin_id])
    test_ompF_required = ompF_array & bin_array   # to see if it has the 3 required mutations
    if sum(test_ompF_required) == 3:
        extra_ompF_required = ompF_array ^ bin_array  # XOR, to check if there are any extra mutations
        if sum(extra_ompF_required) > 0:
            z = 1
    return z


####MAIN BITS####

# Inputs
landscape_input_606 = pd.read_csv('606_wo_extinct.csv')
landscape_input_ecc4 = pd.read_csv('ECC4_wo_extinct.csv')
# death_rate_input = pd.read_csv('death_rate.csv')
# residual_matrix_606 = np.array(pd.read_csv('residuals_606.csv'),dtype=float)
# residual_matrix_ecc4 = np.array(pd.read_csv('residuals_ecc4.csv'),dtype=float)
# -4.20862890184228 for 606 and -7.67936751544442 for ecc4
wildtype_fitness_606 = -4.20862890184228
wildtype_fitness_ecc4 = -7.381
# real_growth_rate of wildtype phage measured by counting plaques in a similar 4hr experiment- ln(100) - ln()
wildtype_fitness_real_606 = np.log(100)
wildtype_fitness_real_ecc4 = np.log(1)
# scaling correction to get real growth rates- ln(100) for 606 and 0 for ecc4
fitness_correction_606 = wildtype_fitness_real_606 - wildtype_fitness_606
fitness_correction_ecc4 = wildtype_fitness_real_ecc4 - wildtype_fitness_ecc4

num_loci = 10

# 629411036 for 606 and 10722464 for ECC4
population_size = 6294110360
time_step = 24*2
generations = 20*time_step
mutation_rate = 7.7E-7
number_runs = 30 # 30
# switch_generation = generations/2

interval=10
switch_time_array = np.array(range(0,generations+int(generations/interval),int(generations/interval)))
total_ompF_switch_time = np.zeros([len(switch_time_array)])

# switch_time_array = switch_time_array[5:]

counter=0;
OmpF_ID = np.array([int(i) for i in srange(1024) if check_ompF(i,10)==1 ])

for switch_generation in switch_time_array:
    
    start=timer()
    
    final_abundance = np.zeros([number_runs, 2**num_loci]);
#     abundance_per_number_mutations = np.zeros([number_runs, num_loci+1]);
#     fitted_selection_rates = np.zeros([2**num_loci,number_runs]);
    population_size_fluctuation = np.zeros([int(generations/time_step),number_runs]);
    population_size_fluctuation[0,:]=population_size 

    print('Total generations {}'.format(generations))
    print('Switch at {}'.format(switch_generation))
    for n in range(0,number_runs,1):
        lambda_landscape_606 = load_landscape(landscape_input_606, wildtype_fitness_606, fitness_correction_606, nearest_neighbour_imputation)#,residual_matrix_606) # this creates an igraph
        filename = 'fittedSelectionRates_new_606_{}.csv'.format(switch_generation)
        save_landscape(lambda_landscape_606, n, filename)
        lambda_landscape_ecc4 = load_landscape(landscape_input_ecc4, wildtype_fitness_ecc4, fitness_correction_ecc4, nearest_neighbour_imputation)#, residual_matrix_ecc4) # this creates an igraph
        filename = 'fittedSelectionRates_new_ECC4_{}.csv'.format(switch_generation)
        save_landscape(lambda_landscape_ecc4, n, filename)
        
        OmpF_abundance = np.zeros([generations, len(OmpF_ID)])
        
        # adds another landscape        
#         lambda_landscape_606_2 = load_landscape(landscape_input_606, wildtype_fitness_606, fitness_correction_606, nearest_neighbour_imputation)#,residual_matrix_606) # this creates an igraph
#         filename = 'fittedSelectionRates_new_606_2_{}.csv'.format(switch_generation)
#         save_landscape(lambda_landscape_606_2, n, filename)
#         lambda_landscape_ecc4_2 = load_landscape(landscape_input_ecc4, wildtype_fitness_ecc4, fitness_correction_ecc4, nearest_neighbour_imputation)#, residual_matrix_ecc4) # this creates an igraph
#         filename = 'fittedSelectionRates_new_ECC4_2_{}.csv'.format(switch_generation)
#         save_landscape(lambda_landscape_ecc4_2, n, filename) 
        
        
### Landscape before switching
        lambda_landscape_1 = lambda_landscape_606   
        lambda_landscape_1['population_size'] = population_size
        lambda_landscape_1.vs[0]['abundance'] = population_size
        
###  Landscape after switching
#         lambda_landscape_2 = lambda_landscape_606_2
        lambda_landscape_2 = lambda_landscape_ecc4
        lambda_landscape_2['population_size'] = population_size

#         plot_population(lambda_landscape_606, 'phage_{0}.png'.format(0))
        print_n = n%10
        [print('Run number {}'.format(n)) if print_n==0 else print_n]
        time_at_switch = 0
#         count_loop = 0
        for t in range(0, switch_generation, time_step): 
            # Change lambda_landscape_xx to ecc4 or 606 in the line below depending on the run
            OmpF_abundance[t:t+time_step,:] = do_timestep(lambda_landscape_1, mutation_rate, time_step, OmpF_ID)#, plot_population, t)
            population_size_fluctuation[int(t/time_step),n] = sum(lambda_landscape_1.vs['abundance'])
            time_at_switch = t
        
        # Change lambda_landscape_xx to ecc4 or 606 in the line below depending on the run
        lambda_landscape_2.vs['abundance'] = lambda_landscape_1.vs['abundance']
#         print('time_at_switch {}'.format(time_at_switch))
#         print('switch_generation {}'.format(switch_generation))
        
        for t in range(time_at_switch+time_step, generations, time_step):
            OmpF_abundance[t:t+time_step,:] = do_timestep(lambda_landscape_2, mutation_rate, time_step, OmpF_ID) 
            population_size_fluctuation[int(t/time_step),n] = sum(lambda_landscape_2.vs['abundance'])      

        final_abundance[n,:] = np.array(lambda_landscape_2.vs['abundance'])    

#         for node in lambda_landscape.vs:
#             abundance_per_number_mutations[n, node['distance']] = abundance_per_number_mutations[n, node['distance']] + node['abundance']


# Saving abundance of all OmpF+ genotypes at all times in a simulation run 
        df = pd.DataFrame(OmpF_abundance)
        filename = 'OmpF-abundance_Switch-generation_{}_run-number_{}.csv'.format(switch_generation,n)
        df.to_csv(filename)
        
    df = pd.DataFrame(population_size_fluctuation)
    filename = 'population-size-fluctuation_Switch-generation_{}_run-number_{}.csv'.format(switch_generation,n)
    df.to_csv(filename)

    df = pd.DataFrame(final_abundance)
    filename = 'final_abundance_{}.csv'.format(switch_generation)
    
    end = timer()
    
    print('Total time for {} generations in this switching - {}'.format(number_runs, end-start) , ' s')
    print('Run time per simulation of switching at this time -',(end-start)/number_runs,'s')

