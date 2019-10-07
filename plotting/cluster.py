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
    energy = np.array(data['movies']['energy'])
    de = energy[1] - energy[0]
    outer_energy = np.linspace(energy[0]-0.5*de, energy[-1]*0.5*de, len(energy)+1)

    n_radial = data['system']['n_radial']
    max_radius = data['system']['max_radius']
    from_center = []
    from_cm = []
    for x in data['collected']:
        from_center.append(x['from_center'])
        from_cm.append(x['from_cm'])
    from_center = np.array(from_center)*1.0;
    from_cm = np.array(from_cm)*1.0;

    r = np.linspace(0, max_radius, n_radial+1)
    R,E = np.meshgrid(r,outer_energy)

    rr = np.linspace(0, max_radius, n_radial)
    for i in range(len(rr)):
        rr[i] = 0.5*(r[i] + r[i+1])
    RR,EE = np.meshgrid(rr,outer_energy)

    # radial /= 4*np.pi*RR**2
    for i in range(from_cm.shape[1]):
        from_center[:,i] /= 4*np.pi/3*(r[i+1]**3 - r[i]**3)
        from_cm[:,i] /= 4*np.pi/3*(r[i+1]**3 - r[i]**3)
    for i in range(from_cm.shape[0]):
        from_center[i,:] *= 1.0/from_center[i,:].max()
        from_cm[i,:] *= 1.0/from_cm[i,:].max()

    plt.figure()
    plt.title(fname)
    plt.pcolor(R, E, from_center)
    plt.colorbar()
    plt.xlabel('$r$')
    plt.ylabel('$E$')

    plt.figure()
    plt.title(fname)
    plt.pcolor(R, E, from_cm)
    plt.colorbar()
    plt.xlabel('$r$')
    plt.ylabel('$E$')

plt.show()
