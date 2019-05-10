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

for my_entropy in sorted(glob.iglob("%s.movie/S*.dat" % filename)):
    my_entropy = my_entropy.replace(' ',',')

dataS = np.loadtxt(my_entropy, ndmin=2)
entropy = dataS[1:,:]
Es = dataS[0,:]
nNs = len(entropy[:,0])
nEs = len(Es)
E_edges = np.zeros(len(Es)+1)
dE = 1

lnw = np.array(yaml_data['bins']['lnw'])
print(np.shape(lnw))
print(lnw , "test")

numE = data['bins']['num_E']
print(numE , "num e")

maxN = data['bins']['max_N']
print(maxN , "MAX N")

EE, NN = np.meshgrid(Es, np.arange(0, maxN+1))

g_exc = np.zeros((maxN+1,numE))
print(EE.shape)
assert(len(lnw) == (maxN+1)*numE)


#wrap the multiplicity from the yaml file into something 2d
for n in range(maxN+1):
    for i in range(0,numE):
        g_exc[n,i] = lnw[n + i*(maxN +1)]
print(g_exc , "mulitplicity")
g_exc -= g_exc[0,0] # set entropy to be zero when there are no atoms.  There is only one such microstate.

print((maxN), 'len of maxN')
print((numE), 'len of numE')
print((nNs), 'len of nNs')
print((nEs), 'len of nEs')


TTT, mumumu = np.meshgrid(np.linspace(0, 10, 50), np.linspace(-19, 9, 50))
NNN = np.zeros_like(TTT)
UUU = np.zeros_like(TTT)
PPP = np.zeros_like(TTT)
TrotterP = np.zeros_like(TTT)
for i in range(TTT.shape[0]):
    for j in range(TTT.shape[1]):
        T = TTT[i,j]
        mu = mumumu[i,j]
        beta = 1/T
        Fid =  NN*T*np.log(NN/V*T**1.5) - NN*T
        Fid[NN==0] = 0
        gibbs_exponent = -beta*(Fid + EE - mu*NN)

        Zgrand = (g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()
        NNN[i,j] = (NN*g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()/Zgrand
        UUU[i,j] = (EE*g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()/Zgrand
        PPP[i,j] = NNN[i,j]*T/V + T/V*(np.log(Zgrand) + np.log(gibbs_exponent.max()*len(gibbs_exponent)))
        TrotterP[i,j] = NNN[i,j]*T/V + T*np.log(Zgrand) / V
    print('mu = {}, T = {}'.format(mu,T))

plt.contourf(TTT, mumumu, NNN, 100)
plt.colorbar()
plt.title('N')
plt.xlabel('T')
plt.ylabel(r'$\mu$')

plt.figure()
plt.contourf(TTT, mumumu, UUU, 100)
plt.colorbar()
plt.title('U')
plt.xlabel('T')
plt.ylabel(r'$\mu$')

plt.figure()
plt.contourf(TTT, mumumu, PPP, 100)
plt.colorbar()
plt.contour(TTT, mumumu, PPP, linewidth=2, color='white')
plt.title('p')
plt.xlabel('T')
plt.ylabel(r'$\mu$')

plt.figure()
plt.contourf(TTT, mumumu, TrotterP, 100)
plt.colorbar()
plt.contour(TTT, mumumu, TrotterP, linewidth=2, color='white')
plt.title('Trotter p')
plt.xlabel('T')
plt.ylabel(r'$\mu$')

TTT, ppp = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 0.00001, 10))
NNN = np.zeros_like(TTT)
UUU = np.zeros_like(TTT)
for i in range(TTT.shape[0]):
    for j in range(TTT.shape[1]):
        T = TTT[i,j]
        p = ppp[i,j]
        beta = 1/T
        Fid =  NN*T*np.log(NN/V*T**1.5) - NN*T
        Fid[NN==0] = 0

        mulo = -10000
        muhi = 10000
        while muhi - mulo > 1e-3:
            mu = 0.5*(muhi + mulo)
            gibbs_exponent = -beta*(Fid + EE - mu*NN)
            Zgrand = (g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()
            pguess_excess = T/V*(np.log(Zgrand) + np.log(gibbs_exponent.max()*len(gibbs_exponent)))
            N_guess = (NN*g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()/Zgrand
            pguess_ideal = N_guess*T/V
            if pguess_excess + pguess_ideal > p:
                muhi = mu
            else:
                mulo = mu
        print('pguess =', pguess_ideal + pguess_excess, pguess_ideal, pguess_excess)

        NNN[i,j] = (NN*g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()/Zgrand
        UUU[i,j] = (EE*g_exc*np.exp(gibbs_exponent - gibbs_exponent.max())).sum()/Zgrand
        print('p = {}, T = {}, mu = {}'.format(p, T, mu))

plt.figure()
plt.contourf(TTT, ppp, UUU, 100)
plt.colorbar()
plt.title('U')
plt.xlabel('T')
plt.ylabel(r'p')

plt.figure()
plt.contourf(TTT, ppp, NNN, 100)
plt.colorbar()
plt.title('N')
plt.xlabel('T')
plt.ylabel(r'p')
plt.show()












