#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from io import StringIO

plt.ion()
for my_histogram in sorted(glob.iglob("samc-1e6-64-movie.yaml/h*.dat")):
    # my_entropy = my_histogram ... switch h to S
    plt.clf()
    print(my_histogram)
    datah = np.loadtxt(my_histogram, ndmin=2)
    print(datah)
    histogram = datah[1:,:]
    Eh = datah[0,:]
    nE = len(Eh)
    nN = len(histogram[:,0])

    N = np.arange(0 , nN+1 , 1)

    plt.pcolor(Eh , N , histogram)
    plt.title('Energy Number Histogram')
    plt.xlabel('Energy')
    plt.ylabel('Number of Atoms')
    plt.colorbar()
    #~ plt.axis([-xL , 0 , 0 , yL])
    plt.pause(.01)
plt.ioff()


plt.figure()
plt.ion()
for my_entropy in sorted(glob.iglob("samc-1e6-64-movie.yaml/S*.dat")):
    plt.clf()
    print(my_entropy)
    dataS = np.loadtxt(my_entropy, ndmin=2)
    print(dataS)
    entropy = dataS[1:,:]
    Es = dataS[0,:]
    nNs = len(entropy[:,0])
    nEs = len(Es)
    #~ temp = np.zeros((nNs,Es))
    #~ for i in range(nNs):
        #~ for j in range(nEs):
            #~ temp[i,j]=(1 / ((entropy[i+1, j] - entropy[i-1,j]) / 4)) 
    #~ print(temp)
    Ns = np.arange(0 , nNs+1 , 1)
    
    E_edges = np.zeros(len(Es)+1)
    dE = 1
    if len(Es) > 1:
        dE = Es[1]-Es[0]
    E_edges[:-1] = Es - 0.5*dE
    E_edges[-1] = Es[-1]
    N_edges = Ns + 0.5
    N_edges[0] = 0
    plt.pcolor(E_edges , N_edges, entropy - entropy[0,-1])
    plt.title('Energy Number Entropy')
    plt.xlabel('Energy')
    plt.ylabel('Number of Atoms')
    plt.colorbar()
    plt.pause(.1)
plt.ioff()

plt.show()

