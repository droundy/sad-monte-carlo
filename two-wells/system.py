import numpy as np
import scipy.special as spl
import sys
import scipy

systems = {
    'easy': {
        'h_small': 1.1,
        'R_small': 0.75,
        'n': 30,
        'min_T': 0.001,
        'min_E': -1.08
    },
    'tiny': {
        'h_small': 1.01,
        'R_small': 0.1,
        'n': 9,
        'min_T': 0.00001,
        'min_E': -1.0099
    },
    'easier': {
        'h_small': 1.1347,
        'R_small': 0.5,
        'n': 9,
        'min_T': 0.001, #FIXME, Figure out what this should be
        'min_E': -1.08 #FIXME, Figure out what this should be
    },
    'easiest': {
        'h_small': 1.0,
        'R_small': 0.5,
        'n': 9,
        'min_T': 0.001, #FIXME, Figure out what this should be
        'min_E': -0.9 #FIXME, Figure out what this should be
    },
    'lj31-like': {
        'h_small': 1.005,
        'R_small': 0.75,
        'n': 90,
        'min_T': 0.001, #FIXME, Figure out what this should be
        'min_E': -1.08 #FIXME, Figure out what this should be
    },
    'hard': {
        'h_small': 1.1,
        'R_small': 0.75,
        'n': 60,
        'min_T': 0.001,
        'min_E': -1.08
    },
    'T-trans-1': {
        'h_small': 1.1,
        'R_small': 0.5,
        'n': 12,
        'min_T': 0.001,
        'min_E': -1.0955
    }
}

if len(sys.argv) == 2:
    system = sys.argv[1]
else:
    system = 'T-trans-1'

min_T = systems[system]['min_T']
min_E = systems[system]['min_E']
h_small = systems[system]['h_small']
R_small = systems[system]['R_small']
n = systems[system]['n']
h_big = 1
R_big = 1.

def V(n):
    return np.pi**(n/2)/(spl.gamma(n/2+1))

def name():
    return system

x_of_cylinder = np.sqrt(R_big**2 - R_small**2)
total_volume = (0.5*V(n)*R_small**n # the small hemisphere
#  + x_of_cylinder*scipy.special.hyp2f1(0.5, (n-1)/2, 1.5, x_of_cylinder**2)
#  - (-R_big)*scipy.special.hyp2f1(0.5, (n-1)/2, 1.5, R_big**2)
  + V(n)*R_big**n # FIXME bad approximation of big sphere with top cut off
)
print('total_volume', total_volume)
# print('hyper', scipy.special.hyp2f1(0.5, (n-1)/2, 1.5, R_big**2))

def D(e):
    if e < -h_big:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2)/total_volume
    return (V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2) + V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2))/total_volume
D = np.vectorize(D)

def D_simplified(e, energy_barrier=0):
    if e < -h_big:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2)/total_volume
    elif e < -h_big+energy_barrier:
        return (V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2) + V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2))/total_volume
    return V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2)

def S(e):
    return np.log(D(e))

def S_simplified(e, energy_barrier=0):
    return np.log(D_simplified(e), energy_barrier=0)
