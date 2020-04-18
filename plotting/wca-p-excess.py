#!usr/bin/python3

import sys, cbor, argparse, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scipy

#Help in Running:
    #currently supports .cbor only
    #   $ python <path to this script> -h
    
#currently, density is read and overwritten as each file is read

parser = argparse.ArgumentParser(description="Create graph of energy vs pressure")
parser.add_argument('cbor', nargs='*',
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


def read_energy(data_loaded, path):
    global my_energy
    my_energy[path] = np.array(data_loaded['movies']['energy'])
    return

def read_entropy(data_loaded, path):
#marks the entropy at each energy at different time t
#movie via entropy shifting for each t; p_exc_tot & energies constant
    global my_entropy
    my_entropy[path] = np.array(data_loaded['movies']['entropy'])
    #for i in range(0, len(my_entropy[path])):
    my_entropy[path][my_entropy==0] = np.nan

def read_pexc_tot(data_loaded, path):
    global my_pexc_tot
    collected = data_loaded['collected']
    my_pexc_tot[path] = np.array([c['pexc_tot'] for c in collected])

def read_count(data_loaded, path):
    global my_count
    collected = data_loaded['collected']
    my_count[path] = np.array([c['count'] for c in collected])

def read_my_t(data_loaded, path):
    global my_t
    my_t[path] = data_loaded['movies']['time']
    for t in range(0, len(my_t[path])):
        my_t[path][t] = float(my_t[path][t])

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
            data_loaded = cbor.load(stream)
        except IOError:
            print('An error occurred trying to read the file.')
    read_energy(data_loaded, path)
    read_my_t(data_loaded, path)
    read_entropy(data_loaded, path)
    read_count(data_loaded, path)
    read_pexc_tot(data_loaded, path)
    read_density(data_loaded)

#i and t are within bounds
def my_temp(path, t, i): #i is the index, t is the time, path the filename
    # dU = TdS - pdV  --->  T = (dU/dS)
    #given E = H - TS, dE/dS=-T hence E/S gives instataneous temp at that condition
    #hence T = - ( (curr-prev E)/(curr - prev S) ) or curr and next E and S
    if i == 0:
        dE = my_energy[path][i+1] - my_energy[path][i]
        dS = my_entropy[path][t][i+1] - my_entropy[path][t][i]
    elif i == len(my_energy[path])-1:
        dE = my_energy[path][i] - my_energy[path][i-1]
        dS = my_entropy[path][t][i] - my_entropy[path][t][i-1]
    else:
        dE = my_energy[path][i+1] - my_energy[path][i-1]
        dS = my_entropy[path][t][i+1] - my_entropy[path][t][i-1]
    try:
        return float(dE/dS)
    except ZeroDivisionError:
        if dE > 0:
            return -float('inf')
        return float('inf') #if dE <= 0

def calc_ideal_p(path, t, i):
    #given pV = NkT where k is boltzmann constant, density is constant so 
    #p = dkt where d is the particular density under examination
#    T = my_temp(t,i)
#    if T!=float('inf') and T!=-float('inf'):
#        return density * scipy.k * T
#    return T #whichever inf it was; +ve or -ve\
    T = my_temperature[path][t][i]
    if T!=float('inf') and T!=-float('inf'):
        return density * scipy.k * T
    return T #whichever inf it was; +ve or -ve

my_pressure={} #my_pressure[path][t][i] -> system's pressure at time t, energy index i of file path
t_size=[] #this array ranks the size of each file's time, as some extend beyond others

for path in args.cbor:
    read_data(path)
    my_temperature[path] = np.zeros_like(my_entropy[path])
    my_pressure[path] = np.zeros_like(my_entropy[path])
    
    #given P = Pexcess + Pideal
    for t in range(0, len(my_t[path])): #the time
        for i in range(0, len(my_energy[path])):  #the energy index
            my_temperature[path][t][i] = my_temp(path,t,i)
            try:
                my_pressure[path][t][i] = my_pexc_tot[path][i] / my_count[path][i] #p_excess
            except ZeroDivisionError:
                if my_pexc_tot[path][i] < 0:
                    my_pressure[path][t][i] = -float('inf')
                else:
                    my_pressure[path][t][i] = float('inf')
                continue
            p_ideal = calc_ideal_p(path, t, i)
            if p_ideal!=float('inf') and p_ideal!=-float('inf'):
                my_pressure[path][t][i] += p_ideal #p_ideal
            else:
                my_pressure[path][t][i] = p_ideal #whatever inf it is; +ve or -ve
   
    #sort t_size - quicksort algorithm as soon as the cbor's data is read
    t_size.append(path)
    s = len(t_size) - 1 # largest index in t_size
    if s >= 1:
        while s > 0:
            curr_path = t_size[s]
            prev_path = t_size[s-1] #the 2nd last path in t_size
            if len( my_t[curr_path] ) < len( my_t[prev_path] ): #if size of my_t of current path greater than previous
                t_size[s] = t_size[s-1]
                t_size[s-1] = curr_path
            else:
                break  #so that worst case is not the average case.                  
        
#remember that t_size[x] gives the path of the cbor, where these are arranged in asending order
#acccording to the length of time they have been running
for t in range(0, len(my_t[ t_size[-1] ])): #uses the timeframe of cbor that has run longest
    plt.figure('energy_temperature')
    plt.clf() #so last figure isn't cleared
    for path in t_size:
        plt.xlabel('Energy (E)')
        plt.ylabel('Temperature (epsilon)')
        plt.title('t = ' + str(my_t[t_size[-1]][t])) #see 2 comments above
        try:
            plt.plot(my_energy[path], my_temperature[path][t], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
        except IndexError:
            plt.plot(my_energy[path], my_temperature[path][-1], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
    plt.pause(0.6)
        
    plt.figure('energy_pressure')
    plt.clf()
    for path in t_size:
        plt.xlabel('Energy (E)')
        plt.ylabel('Pressure (P)')
        plt.title('t = ' + str(my_t[t_size[-1]][t])) #see comments above
        try:
            plt.plot(my_energy[path], my_pressure[path][t], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
        except IndexError:
            plt.plot(my_energy[path], my_pressure[path][-1], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
    plt.pause(0.6)
    
    plt.figure('temperature_pressure')
    plt.clf()
    for path in t_size:
        plt.xlabel('Temperature (epsilon)')
        plt.ylabel('Pressure (P)')
        plt.title('t = ' + str(my_t[t_size[-1]][t])) #see comments above
        try:
            plt.plot(my_temperature[path][t], my_pressure[path][t], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
        except IndexError:
            plt.plot(my_temperature[path][-1], my_pressure[path][-1], label=path)
            plt.tight_layout()
            plt.legend(loc='upper right')
    plt.pause(0.6)
    
plt.show()
