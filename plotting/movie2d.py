#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import math
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from io import StringIO



filename = sys.argv[1]
f = '%s.yaml' % (filename)


# Read YAML file
with open(f, 'r') as stream:
    yaml_data = yaml.load(stream)

data = yaml_data
V = data['system']['cell']['box_diagonal']['x']*data['system']['cell']['box_diagonal']['y']*data['system']['cell']['box_diagonal']['z']

os.system("rm -f tmp*.png") # clean up any preexisting png files

plt.figure()
plt.ion()
i = 0
for my_histogram in sorted(glob.iglob("%s.movie/S*.dat" % filename)):
    # ~ my_entropy = my_histogram ... switch h to S
    plt.clf()
    #print(my_histogram)
    datah = np.loadtxt(my_histogram, ndmin=2)
    #print(datah)
    histogram = datah[1:,:]
    Eh = datah[0,:]
    nE = len(Eh)
    nN = len(histogram[:,0])

    N = np.arange(0 , nN+1 , 1)

    plt.pcolor(Eh , N , histogram)
    plt.xlim(-252,0)
    plt.ylim(0,37)
    plt.title('Energy Number Histogram')
    plt.xlabel('Energy')
    plt.ylabel('Number of Atoms')
    plt.colorbar()
    # ~ # plt.axis([-xL , 0 , 0 , yL])
    plt.pause(.01)
    fname = 'tmp%06d.png'%i # create a file name for this frame
    plt.savefig(fname)
    i += 1

os.system("convert  -delay 50 tmp*.png -loop 2 energy-number-histogram.gif") # make the movie
os.system("rm tmp*.png") # clean up all these png files we created
plt.savefig('energy-number-histogram.png', transparent=True)

plt.figure()
plt.ion()
for my_entropy in sorted(glob.iglob("%s.movie/S*.dat" % filename)):
    my_entropy = my_entropy.replace(' ',',')
    plt.clf()
    #print(my_entropy)
    dataS = np.loadtxt(my_entropy, ndmin=2)
    #print(dataS)
    entropy = dataS[1:,:]
    Es = dataS[0,:]
    print(entropy[0:3,:])
    print(Es)
    nNs = len(entropy[:,0])
    nEs = len(Es)
    # ~ #print("nEs", nEs)
    # ~ #~ #print("nNs",nNs)
    Ns = np.arange(0 , nNs+1 , 1)
    
    E_edges = np.zeros(len(Es)+1)
    dE = 1
    if len(Es) > 1:
        dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
    min_entropy = min(entropy[entropy != 0])
    print('min_entropy', min_entropy)
    if len(entropy) > 1:
        plt.pcolor(E_edges , N_edges, entropy - entropy[1,-1], vmin=-800)
        plt.title('Energy Number Entropy')
        plt.xlabel('Energy')
        plt.ylabel('Number of Atoms')
        plt.colorbar()
        plt.pause(.1)
plt.ioff()  

dataS = np.loadtxt(my_entropy, ndmin=2)
entropy = dataS[1:,:]
Es = dataS[0,:]
nNs = len(entropy[:,0])
nEs = len(Es)
E_edges = np.zeros(len(Es)+1)
dE = 1

plt.figure()
#~ #np.loadtxt(my_entropy[-1], ndim = 2)
temp = np.zeros((nNs,nEs))
for i in range(nNs):
    for j in range(1,nEs-1):
        temp[i,j]=(2 / (entropy[i, j+1] - entropy[i,j-1]))
        #print('temp' , temp[i,j])
        #print('entropy' , entropy[i,j])
print(temp.sum())
if len(Es) > 1:
    dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
print(temp[temp>0])
plt.pcolor(E_edges , N_edges , temp , norm = LogNorm(vmin=.01 , vmax = 100))
plt.title('Temperature')
plt.xlabel('Energy')
plt.ylabel('Number of Atoms')
plt.colorbar()

plt.figure()
almostchempot = np.zeros((nNs,nEs))
chempot = np.zeros((nNs,nEs))
for i in range(1,nNs-1):
    for j in range(nEs):
        #temp[i,j]=(2 / (entropy[i, j+1] - entropy[i,j-1]))
        almostchempot[i,j] = (entropy[i+1,j]-entropy[i-1,j]) / 2
for i in range(1,nNs-1):
    for j in range(1,nEs-1):
        chempot[i,j] = almostchempot[i,j+1] * temp[i+1,j] 
if len(Es) > 1:
    dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
# ~ #print(temp[temp>0])
plt.pcolor(E_edges , N_edges , chempot , norm = LogNorm(vmin=.01 , vmax = 100))
plt.title('Chemical Potential')
plt.xlabel('Energy')
plt.ylabel('Number of Atoms')
plt.colorbar()

plt.figure()
Pexc = np.zeros((nNs,nEs))
for i in range(1,nNs-1):
    for j in range(1,nEs-1):
        Pexc[i,j] = (temp[i+1,j] * entropy[i,j] + chempot[i,j] - j) / 64
if len(Es) > 1:
    dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
# ~ #print(temp[temp>0])nNs,nEs
plt.pcolor(E_edges , N_edges , Pexc , norm = LogNorm(vmin=.01 , vmax = 100))
plt.title('Pressure Excess')
plt.xlabel('Energy')
plt.ylabel('Number of Atoms')
plt.colorbar()



plt.figure()

print('starting meshgrid')
Tmax = 0.5
pmin = 2.5
pmax = 7.0
T_grid, p_grid = np.meshgrid(np.linspace(0, Tmax, 501), np.linspace(0, pmax, 201));

E2d, N2d = np.meshgrid(Es, Ns[:-1])

p_total = Pexc + N2d*temp/64

total_size = E2d.shape[0]*E2d.shape[1]
E_grid = griddata((np.reshape(temp, total_size), np.reshape(p_total, total_size)), np.reshape(E2d, total_size),
                  (T_grid, p_grid), method='nearest')
N_grid = griddata((np.reshape(temp, total_size), np.reshape(p_total, total_size)), np.reshape(N2d, total_size),
                  (T_grid, p_grid), method='nearest')
mu_grid = griddata((np.reshape(temp, total_size), np.reshape(p_total, total_size)), np.reshape(chempot, total_size),
                   (T_grid, p_grid), method='nearest')
print('finished griddata')

plt.pcolor(T_grid, p_grid, mu_grid)
print('finished pcolor')

plt.plot(np.reshape(temp, total_size), np.reshape(p_total, total_size), 'w.')
print('finished dots')


plt.pause(.1)
plt.show()

lnw = np.array(yaml_data['bins']['lnw'])
print(np.shape(lnw))
print(lnw , "test")

numE = data['bins']['num_E']
print(numE , "num e")

maxN = data['bins']['max_N']
print(maxN , "MAX N")

EE, NN = np.meshgrid(Es, np.arange(0, maxN+1))

multi = np.zeros((maxN+1,numE))
print(EE.shape)
assert(len(lnw) == (maxN+1)*numE)
print('multi.shape', multi.shape)


#wrap the multiplicity from the yaml file into something 2d
for n in range(maxN+1):
    for i in range(0,numE):
        multi[n,i] = lnw[n + i*(maxN +1)]
print(multi , "mulitplicity")

print((maxN), 'len of maxN')
print((numE), 'len of numE')
print((nNs), 'len of nNs')
print((nEs), 'len of nEs')


T = 0.5
mu = 1.0
beta = 1/T
Fid = NN*T*np.log(NN/V*T**1.5) - NN*T
gibbs_exponent = -beta*(Fid + EE - mu*NN)

# N_for_this_mu = 

# Zexc = (multi*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()
# print('Zexc', Zexc)


# #beta = 1 / Temp

# #The range needs to change or the indexing in the loop needs to be converted
# #to the actual Energy
# # ~ for n in range(maxN-1):
#     # ~ for e in range(numE):
#         # ~ if not np.isfinite(np.exp((temp[n,e])**-1 * (-e + chempot[n,e] * n))):
#             # ~ zexc = zexc
#         # ~ else:
#             # ~ zexc += np.exp((temp[n,e])**-1 * (-e + chempot[n,e] * n))
# # ~ print(zexc,"Z EXCESS")

# grandU = np.zeros((nNs,nEs))
# zexc = 0
# for i in range(nNs):
#     for j in range(nEs):
#         grandU[i,j] += math.log(((multi[i,j] * j * np.exp((1 / temp[i,j]) * (-j + chempot[i,j] * i)))))
#         zexc += multi[i,j]*(1/(temp[i,j]) * (-j + chempot[i,j] * i))
        
# #~ for E in range(len(Energy)):
#     #~ uexc += gexc[E] * E * math.exp(-beta * E / zexc)

