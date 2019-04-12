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
                           'xkcd:lightblue', 'xkcd:puke']))

my_histogram = {}
current_histogram = {}
my_free_energy = {}
my_volume = {}
current_free_energy = {}
current_total_energy = {}
my_temperature = {}
my_time = {}
my_color = {}
max_iter = 0
my_gamma = {}
my_gamma_t = {}
fnames = sys.argv[1:]
for fname in fnames:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    current_histogram[fname] = np.array(data['bins']['histogram'])
    my_temperature[fname] = data['T']
    current_free_energy[fname] = np.array(data['bins']['lnw'])
    my_volume[fname] = float(data['system']['cell']['box_diagonal']['x'])**3
    current_total_energy[fname] = np.array(data['bins']['total_energy'])
    my_volume[fname] = float(data['system']['cell']['box_diagonal']['x'])**3
    my_color[fname] = allcolors.pop()
    my_time[fname] = np.array(data['movies']['time'])
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    my_temperature[fname] = data['T']
    my_free_energy[fname] = np.array(data['movies']['lnw'])
    my_histogram[fname] = np.array(data['movies']['histogram'])
    my_gamma[fname] = np.array(data['movies']['gamma'], dtype=float)
    my_gamma_t[fname] = np.array(data['movies']['gamma_time'])
    if 'Sad' in data['method']:
        minT = data['method']['Sad']['min_T']


plt.ion()

all_figures = set()
keep_going = True
while keep_going:
    #keep_going = False
    for ii in range(max_iter):
        for fig in all_figures:
            fig.clf()
        for fname in fnames:
            if ii < len(my_time[fname]):
                t = my_time[fname][ii]
                j = ii
            else:
                j = -1

            all_figures.add(plt.figure('Excess free energy'))
            plt.plot(my_free_energy[fname][j,:] - my_free_energy[fname][j,0],
                     my_color[fname],
                     label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$F/kT$')
            plt.legend(loc='best')
            #plt.ylim(Smin, 0)

            all_figures.add(plt.figure('Histogram'))
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.plot(my_histogram[fname][j,:], my_color[fname],
                     label=fname)
            #plt.ylim(0)
            #plt.legend(loc='best')

            all_figures.add(plt.figure('Pressure'))
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('Pressure')
            V = my_volume[fname]
            T = my_temperature[fname]
            F = -my_free_energy[fname][j,:]*T
            N = len(F)
            p = np.zeros(N-1)
            p_exc = np.zeros(N-1)
            for i in range(0,N-1):
                u = F[i+1]-F[i] # dN = 1
                p_exc[i] = (-F[i]+u*(i+.5))/V
                p[i] = (-F[i]+u*(i+.5))/V+(i+.5)*T/V
            UN = np.arange(0.5, N-1, 1)
            print(len(UN), len(p))
            plt.plot(UN, p, my_color[fname],
                     label=fname)
            plt.legend(loc='best')
        plt.figure('Histogram')
        plt.ylabel('histogram')
        plt.ylim(0)
        plt.legend(loc='best')

        plt.pause(0.1)

plt.ioff()
plt.show()
