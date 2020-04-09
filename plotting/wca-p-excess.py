#!usr/bin/python3

import sys, argparse, yaml, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scipy

#Help in Running:
    #currently supports .cbor only
    #   $ python <path to this script> -h

parser = argparse.ArgumentParser(description="Create graph of energy vs pressure")
parser.add_argument('cbor',
                    help = 'path to the .cbor file')
args = parser.parse_args()


#plot of energy vs excess pressure hence:
my_energy={} #the energies
my_entropy={} #entropies[t][i] -> entropy at time t, energy index i
my_pexc_tot={} #excess pressure at energy index i
my_count={} #count of each excess pressure - hence at energy index i
my_t={} #the times, 't'
my_temperature={}  #the temperatures
density = float()


def read_energy(data_loaded):
    global my_energy
    my_energy = np.array(data_loaded['movies']['energy'])
    return

def read_entropy(data_loaded):
#marks the entropy at each energy at different time t
#movie via entropy shifting for each t; p_exc_tot & energies constant
    global my_entropy
    my_entropy = np.array(data_loaded['movies']['entropy'])
    my_entropy[my_entropy==0] = np.nan

def read_pexc_tot(data_loaded):
    global my_pexc_tot
    collected = data_loaded['collected']
    my_pexc_tot = np.array([c['pexc_tot'] for c in collected])

def read_count(data_loaded):
    global my_count
    collected = data_loaded['collected']
    my_count = np.array([c['count'] for c in collected])

def read_my_t(data_loaded):
    global my_t
    my_t = data_loaded['movies']['time']
    for i in range(0, len(my_t)):
        my_t[i] = float(my_t[i])

def read_density(data_loaded):
    global density
    box_diagonal = data_loaded['system']['cell']['box_diagonal']
    positions = data_loaded['system']['cell']['positions']

    volume = box_diagonal['x']*box_diagonal['y']*box_diagonal['z']
    N = len(positions)

    density = N/volume
    return

def read_data(path):
    with open(path, 'rb') as stream:
        try:
            import cbor
            data_loaded = cbor.load(stream)
        except IOError:
            print('An error occurred trying to read the file.')
    read_energy(data_loaded)
    read_entropy(data_loaded)
    read_density(data_loaded)
    read_count(data_loaded)
    read_pexc_tot(data_loaded)
    read_my_t(data_loaded)


#i and t are within bounds
def my_temp(t, i): #i is the index, t is the time
    #given E = H - TS, dE/dS=-T hence E/S gives instataneous temp at that condition
    #hence T = - ( (curr-prev E)/(curr - prev S) ) or curr and next E and S
    if i == 0:
        dE = my_energy[i+1] - my_energy[i]
        dS = my_entropy[t][i+1] - my_entropy[t][i]
    elif i == len(my_energy)-1:
        dE = my_energy[i] - my_energy[i-1]
        dS = my_entropy[t][i] - my_entropy[t][i-1]
    else:
        dE = my_energy[i+1] - my_energy[i-1]
        dS = my_entropy[t][i+1] - my_entropy[t][i-1]
    try:
        return float(-dE/dS)
    except ZeroDivisionError:
        if dE > 0:
            return -float('inf')
        return float('inf') #if dE <= 0

def calc_ideal_p(t, i):
    #given pV = NkT where k is boltzmann constant, density is constant so 
    #p = dkt where d is the particular density under examination
#    T = my_temp(t,i)
#    if T!=float('inf') and T!=-float('inf'):
#        return density * scipy.k * T
#    return T #whichever inf it was; +ve or -ve\
    T = my_temperature[t][i]
    if T!=float('inf') and T!=-float('inf'):
        return density * scipy.k * T
    return T #whichever inf it was; +ve or -ve

read_data(args.cbor)
#given P = Pexcess + Pideal
my_pressure={} #my_pressure[t][i] -> system's pressure at time t, energy index i
for t in range(0, len(my_entropy)): #the time
    my_temperature[t] = np.zeros_like(my_energy)
    my_pressure[t] = np.zeros_like(my_energy)
    for i in range(0, len(my_energy)):  #the index
        my_temperature[t][i] = my_temp(t,i)
        try:
            my_pressure[t][i] = my_pexc_tot[i] / my_count[i] #p_excess
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


for i in range(0, len(my_pressure)):
    plt.figure('pressure')
    plt.clf()
    plt.xlabel('Energy (E)') #need proper symbols
    plt.ylabel('Pressure (P)')
    plt.title('t = ' + str(my_t[i]))
    plt.plot(my_energy, my_pressure[i])

    plt.pause(0.6)

plt.show()
