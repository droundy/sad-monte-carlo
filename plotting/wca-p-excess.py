#!usr/bin/python3

import sys, argparse, yaml, re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#Help in Running:
    #currently supports .cbor only and manually entering reduced density
    #on bash, enter this for further assistance
    #   $ python <path to this script> -h


#QUESTIONS:
#correct to assume gradient at last point equal to the previous (second-last gradient)?
#look at assumed calculations
#what if p_excess is -inf and p_ideal is +inf? P = p_excess + p_ideal
#need proper symbols eg epsilon. Is pressure in pascal, atm etc?
#unsure where to read the reduced density from.


parser = argparse.ArgumentParser(description="Create graph of energy vs pressure")
parser.add_argument('density', type=float,
                    help = 'reduced density')
        #******************density currently in use *****************
parser.add_argument('cbor',
                    help = 'path to the .cbor file')
args = parser.parse_args()


#plot of energy vs excess pressure hence:
Kb=1.38066*(10**-23) #boltzmann constant
my_energy={} #the energies
my_entropy={} #entropies[t][i] -> entropy at time t, energy index i
my_pexc_tot={} #excess pressure at energy index i
my_count={} #count of each excess pressure - hence at energy index i
my_pressure={} #my_pressure[t][i] -> system's pressure at time t, energy index i
my_t={} #the times 't'


def read_energy(data_loaded):
    global my_energy
    my_energy = np.array(data_loaded['movies']['energy'])
    return

def read_entropy(data_loaded):
#marks the entropy at each energy at different time t
#movie via entropy shifting for each t; p_exc_tot & energies constant
    global my_entropy
    my_entropy = data_loaded['movies']['entropy']

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
    #read_density(data_loaded)
    read_count(data_loaded)
    read_pexc_tot(data_loaded)
    read_my_t(data_loaded)


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

read_data(args.cbor)
calc_pressure()

plt.ion()
for i in range(0, len(my_pressure)):
    plt.clf()
    plt.xlabel('Energy\t(E)') #need proper symbols
    plt.ylabel('Pressure\t(P)')
    plt.title('t = ' + str(my_t[i]))
    plt.plot(my_energy, np.array( list(my_pressure[i].keys()) ) )
    if i!=len(my_pressure)-1:
        plt.pause(0.6)
    else:
        plt.show(block=True)
