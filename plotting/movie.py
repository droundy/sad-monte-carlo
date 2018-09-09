#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt

def latex_float(x):
    exp = np.log10(x*1.0)
    if abs(exp) > 2:
        x /= 10.0**exp
        if ('%g' % x) == '1':
            return r'10^{%.0f}' % (exp)
        return r'%g\times 10^{%.0f}' % (x, exp)
    else:
        return '%g' % x

allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']))

my_energy = {}
my_histogram = {}
my_entropy = {}
my_time = {}
my_color = {}
max_iter = 0
my_gamma = {}
my_gamma_t = {}
Smin = None
minT = 0.5
for fname in sys.argv[1:]:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    data['bins']['histogram'] = np.array(data['bins']['histogram'])
    data['bins']['lnw'] = np.array(data['bins']['lnw'])
    my_color[fname] = allcolors.pop()
    my_energy[fname] = np.array(data['movies']['energy'])
    my_time[fname] = np.array(data['movies']['time'])
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    my_entropy[fname] = np.array(data['movies']['entropy'])
    my_histogram[fname] = np.array(data['movies']['histogram'])
    my_gamma[fname] = np.array(data['movies']['gamma'])
    my_gamma_t[fname] = np.array(data['movies']['gamma_time'])
    if 'Sad' in data['method']:
        minT = data['method']['Sad']['min_T']
    if Smin is None:
        Ebest = my_energy[fname];
        Sbest = my_entropy[fname][-1,:]
        Smin = Sbest[Sbest!=0].min() - Sbest.max()

EmaxS = Ebest[np.argmax(Sbest)]
EminT = Ebest[np.argmax(Sbest*minT - Ebest)]
ind_minT = np.argwhere(Ebest == EminT)[0][0]
ind_maxS = np.argwhere(Ebest == EmaxS)[0][0]
print(ind_minT, ind_maxS)
Sbest_interesting = Sbest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]

plt.ion()

plt.figure('gamma')
for fname in my_energy.keys():
        plt.loglog(my_gamma_t[fname], my_gamma[fname], color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'$\gamma$')

plt.figure('comparison')
for fname in my_energy.keys():
    errors = np.zeros(len(my_time[fname]))
    ind_minT = np.argwhere(my_energy[fname] == EminT)[0][0]
    ind_maxS = np.argwhere(my_energy[fname] == EmaxS)[0][0]
    for i in range(len(my_time[fname])):
        S_interesting = my_entropy[fname][i,ind_minT:ind_maxS+1]
        e = S_interesting - Sbest_interesting
        e -= e.mean()
        errors[i] = np.sqrt((e**2).mean())
    plt.loglog(my_time[fname], errors, color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'rms entropy error')

all_figures = set()
while True:
    for i in range(max_iter):
        for fig in all_figures:
            fig.clf()
        all_figures.add(plt.figure('Normed entropy'))
        plt.plot(Ebest, Sbest - Sbest.max(), ':', color='#aaaaaa')
        for fname in my_energy.keys():
            if i < len(my_time[fname]):
                t = my_time[fname][i]
                j = i
            else:
                j = -1

            all_figures.add(plt.figure('Entropy'))
            if j > 0:
                plt.plot(my_energy[fname], my_entropy[fname][j-1,:], my_color[fname],
                         alpha=0.2)
            if j == -1:
                plt.plot(my_energy[fname], my_entropy[fname][j,:], my_color[fname],
                         label=fname+' '+latex_float(len(my_entropy[fname])),
                         alpha=0.2)
            else:
                plt.plot(my_energy[fname], my_entropy[fname][j,:], my_color[fname],
                         label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$S$')
            plt.legend(loc='best')

            all_figures.add(plt.figure('Normed entropy'))
            if j > 0:
                plt.plot(my_energy[fname],
                         my_entropy[fname][i-1,:]-my_entropy[fname][j-1,:].max(),
                         my_color[fname],
                         alpha=0.2)
            if j == -1:
                plt.plot(my_energy[fname],
                         my_entropy[fname][j,:]-my_entropy[fname][j,:].max(),
                         my_color[fname],
                         label=fname+' '+latex_float(len(my_entropy[fname])),
                         alpha=0.2)
            else:
                plt.plot(my_energy[fname],
                         my_entropy[fname][j,:]-my_entropy[fname][j,:].max(),
                         my_color[fname],
                         label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$S$')
            plt.legend(loc='best')
            plt.ylim(Smin, 0)

            all_figures.add(plt.figure('Histogram'))
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('histogram')
            if j > 0:
                plt.plot(my_energy[fname], my_histogram[fname][j-1,:], my_color[fname],
                         alpha=0.2)
            if j == -1:
                plt.plot(my_energy[fname], my_histogram[fname][j,:], my_color[fname],
                         label=fname+' '+latex_float(len(my_entropy[fname])),
                         alpha=0.2)
            else:
                plt.plot(my_energy[fname], my_histogram[fname][j,:], my_color[fname],
                         label=fname)
            plt.legend(loc='best')
        plt.pause(1.0)

plt.ioff()
plt.show()
