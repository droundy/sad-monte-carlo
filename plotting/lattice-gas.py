import yaml
#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

yaml.warnings({'YAMLLoadWarning': False})
#changed from safe_load to full_load to access larger files
with open('ten_by_ten.yaml','r') as stream:
    try:
        data_loaded = yaml.full_load(stream)
    except IOError:
        print('An error occurred trying to read the file.')


#l = 5, E = -50, groups of x in histogram
#L = 100, E = -138, groups of 425 in histogram
#L = 10, E = -34, groups of 13 in histogram

#verified energy for 5x5
#10x10 takes a long time to load, data only appears in last few slides


#template
"""
time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
hist_data = data_loaded['movies']['histogram']
hist = data_loaded['bins']['histogram']
moves = data_loaded['movies']['time']
#energy_data = data_loaded['movies']['energy']

number_data = np.array(data_loaded['movies']['number'])
energy_data = np.array(data_loaded['movies']['energy'])

for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        energy_row = energy_data[i]
        print('row', energy_row)

energy_col = len(energy_data)/energy_row
print('col', energy_col)

#print('energy_resize size', energy_resize.shape)

nlist = len(energy_data)
print('energy_data.shape', energy_data.shape)
energy_data.resize(9, 5)
number_data.resize(9, 5)
print(energy_data)

#num_E +1, Max_N + 1

E = np.zeros((10, 6))
E[:-1,:-1] = energy_data
N = np.zeros((10, 6))
N[:-1,:-1] = number_data

dE = abs(E[0,0] - E[1,0]) #change in energy
E -= dE/2
print('dE', dE)

E[-1,:] = E[-2,:] + dE
N[-1,:] = N[-2,:]

E[:,-1] = E[:,-2]
N[:,-1] = N[:,-2] + 1

N -= 0.5

for t in range(len(entropy_data)):
    #print('time', moves[t])
    S = np.array(entropy_data[t])
    S.resize(9, 5)
    S0 = S[-1,0] # this is the E=0, N=0 entropy
    S = S - S0
    hist = np.array(hist_data[t])
    hist.resize(9, 5)
    S[hist==0] = np.nan
    plt.figure('entropy')
    plt.clf() #Clear the current figure.
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S) #pcolor([X, Y,] C, **kwargs) https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolor.html
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() #https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colorbar.html
    plt.figure('histogram')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,hist)
    plt.colorbar()
    plt.pause(1) #Pause for interval seconds

plt.show()
"""



time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
hist_data = data_loaded['movies']['histogram']
hist = data_loaded['bins']['histogram']
moves = data_loaded['movies']['time']
#energy_data = data_loaded['movies']['energy']

number_data = np.array(data_loaded['movies']['number'])
energy_data = np.array(data_loaded['movies']['energy'])

print(len(energy_data))
#get energy row and col
energy_row = 0
for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        energy_row = int(i + 1)
        #print('row', energy_row)
        break

energy_col = int(len(energy_data)/energy_row)
#print('col', energy_col)

#int(energy_col)
#int(energy_row)


#print('energy_resize size', energy_resize.shape)

nlist = len(energy_data)
#print('energy_data.shape', energy_data.shape)
energy_data.resize(energy_col, energy_row)
number_data.resize(energy_col, energy_row)
#print(energy_data)

"""
Both reshape and resize change the shape of the numpy array;
the difference is that using resize will affect the original array
while using reshape create a new reshaped instance of the array.
"""


E = np.zeros((energy_col + 1, energy_row + 1))
print(E)
E[:-1,:-1] = energy_data
print(E)
N = np.zeros((energy_col + 1, energy_row + 1))
print(N)
N[:-1,:-1] = number_data
print(N)
#why are the E, N arrays sliced?


#E = energy_data
#N = number_data

dE = abs(E[0,0] - E[1,0]) 
E -= dE/2
print('dE', dE)

E[-1,:] = E[-2,:] + dE
N[-1,:] = N[-2,:]

E[:,-1] = E[:,-2]
N[:,-1] = N[:,-2] + 1

N -= 0.5


for t in range(len(entropy_data)):
    #print('time', moves[t])
    S = np.array(entropy_data[t])
    S.resize(energy_col, energy_row)
    S0 = S[-1,0] # this is the E=0, N=0 entropy
    S = S - S0
    hist = np.array(hist_data[t])
    hist.resize(energy_col, energy_row)
    S[hist==0] = np.nan
    plt.figure('entropy')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S) #pcolor([X, Y,] C, **kwargs) https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolor.html
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() #https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colorbar.html
    plt.figure('histogram')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,hist)
    plt.colorbar()
    plt.pause(1)

plt.show()
