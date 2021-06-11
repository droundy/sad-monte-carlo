# -*- coding: utf-8 -*-
"""
Created on Sun May 23 21:35:27 2021

@author: holme
"""
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
from scipy.special import erf
from scipy.special import erfi
'''
ETROPY OF FLUIDS - XENO'S METHOD

'''

'''
Gaussian Test Data
'''
# define varaiables
dE = 0.01
E = np.arange(-5,5,0.01)
Eb = np.array([5,2,1,0,-1,-2,-3]) # Eborders
sigma = 1
N = len(Eb) # number of borders


# initialize arrays
D = np.zeros((len(E),1))
w = np.zeros((N-1,1))
E_avg = np.zeros((N-1,1))

# calculate density of states
for i in range(len(E)-1):
    D[i] = np.exp(-E[i]**2/(2*sigma**2))/np.sqrt(np.pi)

# calculate weights for each bin
for i in range(1,N):
    w[i-1] = -(erf((Eb[i]/np.sqrt(2))) - erf(Eb[i-1]/(np.sqrt(2))))/2

# calculate average energy for each bin
for i in range(1,N):
    E_avg[i-1] = (np.exp((-E[i]**2)/2)-np.exp((-E[i-1]**2)/2))/(w[i-1]*np.sqrt(2*np.pi))

print('w = ', w)
print('the sum of the weights is ', sum(w))   
print('E = ', E_avg)
print('the sum of the avg Energies is ', sum(E_avg))
    
# # plot

plt.plot(E,D)
# for e in Eb:
#     plt.axvline(e,color='k')
    
plt.ylabel('Density of States')
plt.xlabel('Energy')

x = [-2.5,-1.5,-0.5,0.5,1.5,2.5]
plt.plot(x,w,'o', color = 'blue')
plt.plot(x,w,'--', color = 'blue')
# plt.plot(x,E_avg,'o',color = 'red')
# plt.plot(x,E_avg,'--',color = 'red')
# plt.plot(x,E_avg**2, 'o', color = 'lime')
# plt.plot(x,E_avg**2, '--', color = 'lime')
# plt.xlim(-3,3)
    
# plt.show()


'''
Xeno's Method
'''

### Root Finding algorithm
def Roots(A):
    pass

### Calculate Methods

# Notes:
# SL = Si
# SR = Si-1
# A must be sorted into positive or negative before applying erfi or erf

def Moments(SL, SR, A, EL, ER): 
    # complex variables
    epsilon = ER - EL
    gamma = A*epsilon**2
    delta_S = SR - SL
    
    # calculate methods
    if A < 0:
        weight = (np.sqrt(np.pi)*np.exp(SL)*np.exp((gamma-delta_S)**2/4*gamma)*(erf((delta_S + gamma)/(2*np.sqrt(gamma)))-erf((delta_S - gamma)/(2*np.sqrt(gamma)))))/(2*np.sqrt(gamma))
        avg_E = epsilon*(((delta_S + gamma)/(2*gamma))-((1-np.exp(-delta_S))/(erf((delta_S + gamma)/2*np.sqrt(gamma))-erf((delta_S + gamma)/2*np.sqrt(gamma))))) + EL
        theta = ((np.sqrt(np.pi)/8*gamma**(5/2))((delta_S + gamma)**2+2*gamma))*(-np.exp((delta_S + gamma)**2/4*gamma + delta_S))*erf((delta_S + gamma)/2*np.sqrt(gamma)) + \
            np.exp(((delta_S - gamma)**2)/4*gamma)*erf((delta_S - gamma)/2*np.sqrt(gamma)) + \
            - ((3*gamma - delta_S)-np.exp(delta_S)*(gamma - delta_S))/(4*gamma**2)
        avg_E_2 = (epsilon*np.exp(SL)*theta)/weight + 3*EL**2 - 2*EL*avg_E
        stdev_E = np.sqrt(avg_E**2 - avg_E_2)
        return [weight,avg_E,avg_E_2,stdev_E]
    elif A > 0:
        weight = (np.sqrt(np.pi)*np.exp(SL)*np.exp((gamma-delta_S)**2/4*gamma)*(erfi((delta_S + gamma)/(2*np.sqrt(gamma)))-erfi((delta_S - gamma)/(2*np.sqrt(gamma)))))/(2*np.sqrt(gamma))
        avg_E = epsilon*(((delta_S + gamma)/(2*gamma))-((1-np.exp(-delta_S))/(erfi((delta_S + gamma)/2*np.sqrt(gamma))-erfi((delta_S + gamma)/2*np.sqrt(gamma))))) + EL
        theta = ((np.sqrt(np.pi)/8*gamma**(5/2))((delta_S + gamma)**2+2*gamma))*(-np.exp((delta_S + gamma)**2/4*gamma + delta_S))*erfi((delta_S + gamma)/2*np.sqrt(gamma)) + \
            np.exp(((delta_S - gamma)**2)/4*gamma)*erfi((delta_S - gamma)/2*np.sqrt(gamma)) + \
            - ((3*gamma - delta_S)-np.exp(delta_S)*(gamma - delta_S))/(4*gamma**2)
        avg_E_2 = (epsilon*np.exp(SL)*theta)/weight + 3*EL**2 - 2*EL*avg_E
        stdev_E = np.sqrt(avg_E**2 - avg_E_2)
        return [weight,avg_E,avg_E_2,stdev_E]
    else:
        return "A == 0"

# function I want:
def FindA_SL_and_SR(known_weight, known_avg_E, known_avg_E_2, EL, ER):
    # guess initial SR, SL, A, then use root? to find them accurately

    # fucntion I need for root-finding (or use lambda expression)
    def find_residual(A, SL, SR):
        w, avg_E, avg_E_2, stdev_E = Moments(A, SL, SR, EL, ER)
        return [w-known_weight, avg_E - known_avg_E, etc]

    return A, SL, SR
### Compile into 1 funciton

def Xeno(SL, SR, A, EL, ER):
    pass

### Plot Moments

    

                              