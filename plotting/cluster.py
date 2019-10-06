#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt

def latex_float(x):
    exp = int(np.log10(x*1.0))
    if abs(exp) > 2:
        x /= 10.0**exp
        if ('%.1g' % x) == '1':
            return r'10^{%.0f}' % (exp)
        return r'%.1g\times10^{%.0f}' % (x, exp)
    else:
        return '%g' % x

allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
                           'xkcd:lightblue', 'xkcd:puke', 'xkcd:puce', 'xkcd:turquoise']))

my_energy = {}
my_histogram = {}
my_entropy = {}
my_time = {}
my_color = {}
max_iter = 0
my_gamma = {}
my_gamma_t = {}
Smin = None
minT = 1.0
fnames = sys.argv[1:]
for fname in fnames:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    data['system']['bins']['radial'] = np.array(data['system']['bins']['radial'])
    radial = data['system']['bins']['radial']
    radial = radial*1.0
    r = np.linspace(0, data['system']['max_radius'], data['system']['bins']['n_radial']+1)
    width = data['system']['bins']['width']
    mine = data['system']['bins']['min']
    e = np.arange(mine,
                  mine + 0.5*width + radial.shape[0]*width,
                  width)
    R,E = np.meshgrid(r,e)

    rr = np.linspace(0, data['system']['max_radius'], data['system']['bins']['n_radial'])
    for i in range(len(rr)):
        rr[i] = 0.5*(r[i] + r[i+1])
    ee = np.arange(mine + 0.5*width,
                   mine + radial.shape[0]*width,
                   width)
    RR,EE = np.meshgrid(rr,ee)

    radial /= 4*np.pi*RR**2
    for i in range(radial.shape[0]):
        radial[i,:] *= 1.0/radial[i,:].max()
    plt.pcolor(R, E, radial)
    plt.colorbar()
    plt.xlabel('$r$')
    plt.ylabel('$E$')

plt.show()
