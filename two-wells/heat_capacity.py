#!/usr/bin/python3
import os
import numpy as np
import matplotlib.pyplot as plt
import system, compute
import glob

def C(T, S):#T is a temperature and S is an entropy function
    E = np.linspace(-system.h_small, 0, 5000)
    E = 0.5*(E[1:] + E[:-1])
    dE = E[1] - E[0]

    def normalize_S(S):
        S = S - max(S)
        total = np.sum(np.exp(S)*dE)
        return S - np.log(total)

    S = normalize_S(S(E))

    S_minus_E = S-E/T
    M = np.max(S_minus_E)

    Z = np.sum(np.exp(S_minus_E-M))*dE

    avg_E = np.sum(np.exp(S_minus_E-M) * E) * dE / Z

    avg_E_squared = np.sum(np.exp(S_minus_E-M) * E**2) * dE / Z

    return (avg_E_squared - avg_E**2 ) / T**2

#Testing
#t = np.linspace(0.005,0.012,100)
#c = [C(T) for T in t]
#plt.plot(t, c)
#plt.show()







