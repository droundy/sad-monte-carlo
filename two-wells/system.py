import numpy as np
import scipy.special as spl
import sys

systems = {
    'easy': {
        'h_small': 1.1,
        'R_small': 0.5,
        'n': 30,
    },
    'easier': {
        'h_small': 1.1347,
        'R_small': 0.5,
        'n': 9,
    },
    'lj31-like': {
        'h_small': 1.005,
        'R_small': 0.75,
        'n': 90,
    }
}

if len(sys.argv) == 2:
    system = sys.argv[1]
else:
    system = 'easier'

h_small = systems[system]['h_small']
R_small = systems[system]['R_small']
n = systems[system]['n']
h_big = 1
R_big = 1.

def V(n):
    return np.pi**(n/2)/(spl.gamma(n/2+1))

def D(e):
    if e < -h_big:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2)
    return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2) + V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2)
D = np.vectorize(D)

def S(E):
    return np.log(D(E))