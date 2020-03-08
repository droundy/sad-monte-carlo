#!usr/bin/python3

import sys, argparse, yaml, re
import numpy as np
import matplotlib.pyplot as plt

#gradient at last point equal to the previous (second-last gradient)
#look at assumed calculations
#what is p_excess is -inf and p_ideal is +inf?
#number of entropies and energies unequal. how were they stored?

parser = argparse.ArgumentParser(description="Create graph of energy vs pressure")
parser.add_argument('density', type=float,
                    help = 'reduced density')
parser.add_argument('cbor',
                    help = 'path to the .cbor file')
args = parser.parse_args()



#plot of energy vs excess pressure hence:
Kb=1.38066*(10**-23) #boltzmann constant ###
my_energy={} #array of the energies read directly from .energy file ###
my_entropy={} #array of entropies ###
my_pexc_tot={} #excess pressure  ###
my_count={} #count at each excess pressure ###
my_pressure={} #pressure of system
my_histogram={} #histogram data

def read_energy(data_loaded):
    global my_energy
    my_energy = np.array(data_loaded['movies']['energy'])
    return

def read_entropy(data_loaded):
#x is all 151 - lines up with rest of data
#marks the entropy at each energy at different time t
#movie via entropy shifting for each t; p_exc_tot & energies constant
    global my_entropy
    my_entropy = data_loaded['movies']['entropy']
    return

def read_pexc_tot(data_loaded):
    global my_pexc_tot
    collected = data_loaded['collected']
    my_pexc_tot = np.array([c['pexc_tot'] for c in collected])

def read_count(data_loaded):
    global my_count
    collected = data_loaded['collected']
    my_count = np.array([c['count'] for c in collected])

def read_density(data_loaded):
    return

def read_hist(data_loaded):
    return

def read_data(path):
    with open(path) as stream:
        try:
            import cbor
            data_loaded = cbor.load(stream)
        except IOError:
            print('An error occurred trying to read the file.')
    read_energy(data_loaded)
    read_entropy(data_loaded)
    read_density(data_loaded)
    read_hist(data_loaded)
    read_count(data_loaded)
    read_pexc_tot(data_loaded)
    return


#i and t are within bounds
def my_temp(t, i): #i is the index, t is the time
    #given E = H - TS, dE/dS=-T hence E/S gives instataneous temp at that condition
    #hence T = - ( (curr-prev E)/(curr - prev S) ) or curr and next E and S
    if i == len(my_energy)-1: #this assumes last gradient similar to previous gradient
        i-=1
    dE = float(my_energy[i+1]) - float(my_energy[i])
    dS = float(my_entropy[t][i+1]) - float(my_entropy[t][i])
    try:
        return float(-dE/dS)
    except ZeroDivisionError:
        if dE > 0:
            return -float('inf')
        return float('inf') #if dE <= 0

def calc_ideal_p(t, i):
    #given pV = NkT where k is boltzmann constant, density is constant so 
    #p = dkt where d is the particular density under examination
    T = my_temp(t,i)
    if T!=float('inf') and T!=-float('inf'):
        return args.density * Kb * my_temp(t,i)
    return T #whichever inf it was; +ve or -ve

#function assumes at a given index i, my_energy[i] corresponds with my_entropy[i]
#as well as my_pexc_tot[i] and my_count[i]
def calc_pressure():   #calculate pressure and populate the my_pressure array
    #given P = Pexcess + Pideal
    global my_pressure
    for t in range(0, len(my_entropy)):
        my_pressure[t]={}
        for i in range(0, len(my_energy)):
            try:
                my_pressure[t][i] = float(my_pexc_tot[i]) / float(my_count[i]) #p_excess
            except ZeroDivisionError:
                if my_pexc_tot[i] < 0:
                    my_pressure[t][i] = -float('inf')
                else:
                    my_pressure[t][i] = float('inf')
                continue
            p_ideal = calc_ideal_p(t, i)
            if p_ideal!=float('inf') and p_ideal!=-float('inf'):
                my_pressure[t][i] += calc_ideal_p(t, i) #p_ideal
            else:
                my_pressure[t][i] = p_ideal #whatever inf it is; +ve or -ve
    return


read_data(args.cbor)
calc_pressure()
#plot energy vs pressure

#t=1
#for i in range(0, len(my_energy)): 
#    try:
#        out = float(my_pexc_tot[i]) / float(my_count[i]) #p_excess
#    except ZeroDivisionError:
#        if my_pexc_tot[i] < 0:
#            out = -float('inf')
#        else:
#            out = float('inf')
#        #print(out)
#        continue
#    p_ideal = calc_ideal_p(t, i)
#    if p_ideal!=float('inf') and p_ideal!=-float('inf'):
#        out += calc_ideal_p(t, i) #p_ideal
#    else:
#        out = p_ideal #whatever inf it is; +ve or -ve
#    print(out)

#for i in range(0, len(my_energy)): 
#    try:
#        print(float(my_pexc_tot[i]) / float(my_count[i])) #p_excess
#    except ZeroDivisionError:
#        if my_pexc_tot[i] < 0:
#            print(-float('inf'))
#        print(float('inf'))


#fname = args.yaml
#with open(args.yaml) as y:
#    yaml_data = y.read()
#    data = yaml.load(yaml_data)
#    collected = data['collected']
#    my_pexc_tot = np.array([c['pexc_tot'] for c in collected])
#    my_count = np.array([c['count'] for c in collected])
