#!usr/bin/python3

import sys, argparse, yaml, re
import numpy as np
import matplotlib.pyplot as plt

#what is the boltzmann constant
#look at assumed calculations
#there were '.inf' values in the entropy file. what are these?
#number of entropies and energies unequal. how were they stored?

parser = argparse.ArgumentParser(description="Create graph of energy vs pressure")
parser.add_argument('density', type=float,
                    help = 'reduced density')
parser.add_argument('yaml',
                    help = 'path to the .yaml file')
parser.add_argument('energy',
                    help = 'path to the .energy file')
parser.add_argument('entropy',
                    help = 'path to the .entropy file')
args = parser.parse_args()



#plot of energy vs excess pressure hence:
Kb=1.38066*(10**-23) #boltzmann constant    ******************
my_energy={} #array of the energies read directly from .energy file ###
my_entropy={} #array of entropies
my_pexc_tot={} #excess pressure  ###
my_count={} #count at each excess pressure ###
my_pressure={} #pressure of system
my_histogram={} #histogram data


#The following functions (read*:) dependent on maintainance of current storage method
#of the data in .energy, .entropy etc files.
def read_energy_data(path):
    with open(path) as f:
        temp = (f.read().split("\n"))
    for i in range(0, len(temp)-1): #exclude \n at end
#in this case, the splitting created empty array location due to \n at end. Thus -1 is used
#to avoid that empty location (array bounds not of key consideration).
        my_energy[i] = float(temp[i])
    return

#will use this or histogram to eliminate div by 0
def read_entropy_data(path):
    with open(path) as f:
        temp = f.read().split("\t")
    for i in range(0, len(temp)):
        temp[i] = (temp[i].split("\n")[0])
    for i in range(0, len(temp)):
        my_entropy[i] = float(temp[i])
    return

def my_temp(e, s):
    return
    #given E = H - TS, dE/dS=-T hence E/S gives instataneous temp at that condition
    #hence T = - ( (curr-prev E)/(curr - prev S) ) or curr and next E and S

def calc_ideal_p(e, s):
    #given pV = NkT where k is boltzmann constant, density is constant so 
    #p = dkt where d is the particular density under examination
    return args.density * Kb * my_temp(e,s)

#function assumes at a given index i, my_energy[i] corresponds with my_entropy[i]
#as well as my_pexc_tot[i] and my_count[i]
def calc_p():   #calculate pressure and populate the my_pressure array
    #given P = Pexcess + Pideal
    p_ideal = {}
    for i in range(0, len(my_energy)):
        p_ideal[i] = calc_ideal_p(my_energy[i], my_entropy[i])
        my_pressure[i] = ( my_pexc_tot[i] / my_count[i] ) + p_ideal[i]

fname = args.yaml
with open(args.yaml) as y:
    yaml_data = y.read()
    data = yaml.load(yaml_data)
    collected = data['collected']
    my_pexc_tot = np.array([c['pexc_tot'] for c in collected])
    my_count = np.array([c['count'] for c in collected])
    entropy = np.array(data['movies']['entropy'])
    print(entropy.shape)

read_energy_data(args.energy)
read_entropy_data(args.entropy) #** need to pop .inf if needed

#calc_p()
#plot the pressure vs E
