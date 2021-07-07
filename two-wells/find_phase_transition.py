#!/usr/bin/python3

import numpy as np
import scipy.special as spl
import matplotlib.pyplot as plt

h_small = 1.005
h_big = 1
R_small = 0.75
R_big = 1.
n = 90

E1 = -133.58642 # minimum energy (Mackay) for an LJ31 cluster
E2 = -133.29382 # first local minimum (anti-Mackay) for an LJ31 cluster
E_barrier = -131 # approximate energy of the transition state between the two
T_transition = 0.025 # approximate temperature for transition between the two
delta_E = E2 - E1

print('temperature over Delta E LJ31', T_transition/delta_E)
print('barrier over Delta E LJ31', (E_barrier-E1)/delta_E)


def V(n):
    return np.pi**(n/2)/(spl.gamma(n/2+1))

def D(e):
    if e < -h_big:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2)
    else:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2) + V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2)


def E_1(e):
    '''
       Finds the lower energy that has the same slope as energy e
    '''
    numerator = R_small**n * np.sqrt(e/h_small+1)**(n-2) + R_big**n * np.sqrt(e/h_big+1)**(n-2)

    denominator = (1/h_small)*R_small**n * np.sqrt(e/h_small+1)**(n-4) + (1/h_big)*R_big**n * np.sqrt(e/h_big+1)**(n-4) 

    return numerator/denominator - h_small

def S(e):
    return np.log(D(e)/D(0))

def S_prime(e): #ONLY VALID FOR VALUES BELOW SHORT WELL. More complicated above the short well
    return ((n-2)/(2*h_small))/(e/h_small + 1)

def g(e):
    return S(e) - S(E_1(e)) - S_prime(E_1(e))*(e - E_1(e))

##### Bisection method to find zero of g #####
max_iter = 1000


lower_bound = -h_big#h_2*0.99
upper_bound = -0.1*h_big #The highest E_2 with a valid E_1
guess = (lower_bound + upper_bound)/2 #Guess based on looking at the plot

iters = 0
print(lower_bound, '< E000 <', upper_bound)
while upper_bound - lower_bound > 1e-10 and iters < max_iter:
    iters += 1
    if g(guess) > 0 or np.isnan(g(guess)) or abs(E_1(guess) - guess) < 1e-10:
        upper_bound = guess
    else:
        lower_bound  = guess
    guess = (lower_bound + upper_bound) / 2

print('\nguess is', guess)

# print("Guess: " + str(guess) )
# print("g(guess): " +  str(g(guess)) )
# print("1/T: " + str(S_prime(guess)))
# print("->T: " + str(1/S_prime(guess)))

actual_T = 1/S_prime(guess)
print('actual T/Delta E', actual_T/(h_small - h_big))
print('actual barrier/Delta E', h_small/(h_small - h_big))

print('barrier energy should be', (actual_T/0.025)*(E_barrier-E1) - h_small)

##### Plotting #####
print('guess', guess, E_1(guess))
E_2 = np.arange(-h_big, 0, 0.0001)


# plt.figure()
# plt.plot(E_2, E_1(E_2))
# plt.plot([-h_big,0],[-h_big,-h_big])
# plt.title(r'$E_1(E_2)$')

e_2 = []

for e in E_2:
    if E_1(e) < -h_big:
        e_2.append(e)

e_1 = []
for e in e_2:
    e_1.append(E_1(e))

plt.figure()
E = np.arange(-h_small,0,abs(h_big-h_small)/100)

s = np.zeros_like(E)

for i in range(0,len(E)):
    s[i] = S(E[i])

# s_1 = []
# s_2 = []
# for i in range(0,len(e_1)):
#     s_1.append(np.log( D(e_1[i])) )
#     s_2.append(np.log( D(e_2[i])) )

plt.plot(E,s)
plt.plot(guess, S(guess), '.')
plt.plot(E_1(guess), S(E_1(guess)), '.')
# plt.plot(e_1,s_1)
# plt.plot(e_2,s_2)

##### Common Tangent #####

t = np.arange(-h_small,-0,0.5)
m = S_prime(E_1(guess))
print('T =', 1/m)
Y = m*(t-guess) + S(guess)

plt.plot(t,Y)
plt.ylim(ymax=0)
plt.title(r'$S(E)$')
##### g(E_2) #####

# G = np.zeros_like(E_2) 

# for i in range(0,len(G)):
#     G[i] = g(E_2[i])
# plt.figure()
# plt.plot(E_2,G)


# plt.title(r'$G(E_2)$')
plt.show()