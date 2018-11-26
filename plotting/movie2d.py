#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.colors import LogNorm
from io import StringIO


filename = sys.argv[1]
f = '%s.yaml' % (filename)


# Read YAML file
with open(f, 'r') as stream:
    yaml_data = yaml.load(stream)

data = yaml_data
lnw = data['bins']['lnw'] = np.array(yaml_data['bins']['lnw'])
#print(np.shape(lnw))


#~ plt.ion()
#~ for my_histogram in sorted(glob.iglob("samc-1e6-64-movie.yaml/h*.dat")):
    #~ # my_entropy = my_histogram ... switch h to S
    #~ plt.clf()
    #~ #print(my_histogram)
    #~ datah = np.loadtxt(my_histogram, ndmin=2)
    #~ #print(datah)
    #~ histogram = datah[1:,:]
    #~ Eh = datah[0,:]
    #~ nE = len(Eh)
    #~ nN = len(histogram[:,0])

    #~ N = np.arange(0 , nN+1 , 1)

    #~ plt.pcolor(Eh , N , histogram)
    #~ plt.title('Energy Number Histogram')
    #~ plt.xlabel('Energy')
    #~ plt.ylabel('Number of Atoms')
    #~ plt.colorbar()
  #  plt.axis([-xL , 0 , 0 , yL])
    #~ plt.pause(.01)
#~ plt.ioff()


#~ plt.figure()
#~ plt.ion()
for my_entropy in sorted(glob.iglob("samc-1e6-64-movie.yaml/S*.dat")):
    my_entropy = my_entropy.replace(' ',',')
    #~ #print(my_entropy)
    plt.clf()
    #print(my_entropy)
    dataS = np.loadtxt(my_entropy, ndmin=2)
    #print(dataS)
    entropy = dataS[1:,:]
    Es = dataS[0,:]
    nNs = len(entropy[:,0])
    nEs = len(Es)
    #~ print("nEs", nEs)
    #~ print("nNs",nNs)
    Ns = np.arange(0 , nNs+1 , 1)
    
    E_edges = np.zeros(len(Es)+1)
    dE = 1
    if len(Es) > 1:
        dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
    #~ plt.pcolor(E_edges , N_edges, entropy - entropy[0,-1])
    #~ plt.title('Energy Number Entropy')
    #~ plt.xlabel('Energy')
    #~ plt.ylabel('Number of Atoms')
    #~ plt.colorbar()
    #~ plt.pause(.1)
#~ plt.ioff()  

dataS = np.loadtxt(my_entropy, ndmin=2)
entropy = dataS[1:,:]
Es = dataS[0,:]
nNs = len(entropy[:,0])
nEs = len(Es)
E_edges = np.zeros(len(Es)+1)
dE = 1

#~ np.loadtxt(my_entropy[-1], ndim = 2)
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

E2d, N2d = np.meshgrid(Es, Ns[:-1])
plt.figure()
print(temp.shape, N2d.shape, E2d.shape)
plt.contourf(temp.T, N2d.T, E2d.T)
plt.title('Energy')
plt.xlabel('T')
plt.ylabel('Number of Atoms')
plt.xlim(0, 5)
plt.ylim(0, nNs)
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
#print(temp[temp>0])
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
#print(temp[temp>0])
plt.pcolor(E_edges , N_edges , Pexc , norm = LogNorm(vmin=.01 , vmax = 100))
plt.title('Pressure Excess')
plt.xlabel('Energy')
plt.ylabel('Number of Atoms')
plt.colorbar()

plt.pause(.1)
plt.show()

