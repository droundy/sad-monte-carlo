import numpy as np
import scipy.special as spl
import sys
import scipy

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


x_of_cylinder = np.sqrt(R_big**2 - R_small**2)
total_volume = (0.5*V(n)*R_small**n # the small hemisphere
 + V(n-1)*R_small**(n-1)*(R_small + R_big - x_of_cylinder) # the small cylinder
 + x_of_cylinder*scipy.special.hyp2f1(0.5, (n-1)/2, 1.5, x_of_cylinder**2)
 - (-R_big)*scipy.special.hyp2f1(0.5, (n-1)/2, 1.5, R_big**2)
)

def D(e):
    if e < -h_big:
        return V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2)/total_volume
    return (V(n)*R_small**n/2*np.sqrt(e/h_small+1)**(n-2) + V(n)*R_big**n/2*np.sqrt(e/h_big+1)**(n-2))/total_volume
D = np.vectorize(D)

def S(E):
    return np.log(D(E))