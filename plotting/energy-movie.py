#!/usr/bin/python3

import yaml, sys, argparse
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

parser = argparse.ArgumentParser(description="create movie and graphs for energy histogram data")
parser.add_argument('--minT', action='store', type=float, default=0.005,
                    help = "the minimum temperature of interest")
parser.add_argument('--maxT', action='store', type=float, default=2.0,
                    help = "the maximum temperature of interest")
parser.add_argument('--match-energy', action='store', type=float,
                    help = "the energy at which we want to normalize the entropy")
parser.add_argument('yaml', nargs='*',
                    help = 'the names of some yaml files')
args = parser.parse_args()

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
minT = 1.1
for fname in args.yaml:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    print('Done loading yaml')
    data['bins']['histogram'] = np.array(data['bins']['histogram'])
    data['bins']['lnw'] = np.array(data['bins']['lnw'])
    my_color[fname] = allcolors.pop()
    my_energy[fname] = np.array(data['movies']['energy'])
    my_time[fname] = np.array(data['movies']['time'])
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    my_entropy[fname] = np.array(data['movies']['entropy'])
    if args.match_energy is not None:
        match_index = (np.abs(my_energy[fname] - args.match_energy)).argmin()
        for i in range(len(my_entropy[fname][:,match_index])):
            my_entropy[fname][i,:] -= my_entropy[fname][i,match_index]

    my_histogram[fname] = np.array(data['movies']['histogram'])
    my_gamma[fname] = np.array(data['movies']['gamma'])
    my_gamma_t[fname] = np.array(data['movies']['gamma_time'])
    if 'Sad' in data['method']:
        minT = data['method']['Sad']['min_T']
    if Smin is None:
        Ebest = my_energy[fname];
        Sbest = my_entropy[fname][-1,:]
        hbest = my_histogram[fname][-1,:]
        Smin = Sbest[Sbest!=0].min() - Sbest.max()
        if args.match_energy is not None:
            Smin = Sbest[hbest!=0].min()

ind_maxS = np.argmax(Sbest)
EmaxS = Ebest[ind_maxS]
ind_minT = np.argmax(Sbest*minT - Ebest)
EminT = Ebest[ind_minT]
print('energies:', Ebest[ind_minT], 'at temperature', minT, 'and max entropy', Ebest[ind_maxS])
Sbest_interesting = Sbest[ind_minT:ind_maxS+1]
Ebest_interesting = Ebest[ind_minT:ind_maxS+1]

plt.ion()

plt.figure('gamma')
for fname in args.yaml:
        plt.loglog(my_gamma_t[fname], my_gamma[fname], color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'$\gamma$')
# plt.ylim(1e-12, 1.1)

def heat_capacity(T, E, S):
    C = np.zeros_like(T)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P = P/P.sum()
        U = (E*P).sum()
        C[i] = ((E-U)**2*P).sum()/T[i]**2
    return C

plt.figure('comparison')
for fname in args.yaml:
    if my_energy[fname][0] <= EminT:
        errors = np.zeros(len(my_time[fname]))
        ind_minT = np.argmin(np.abs(my_energy[fname] - EminT))
        ind_maxS = np.argmin(np.abs(my_energy[fname] - EmaxS))
        for i in range(len(my_time[fname])):
            # We interpolate to find the values of the entropy at the
            # energies that show up in our "best" data.
            S_interesting = np.interp(Ebest_interesting, my_energy[fname], my_entropy[fname][i,:])
            e = S_interesting - Sbest_interesting
            e -= e[Sbest_interesting!=0.0].mean()
            errors[i] = np.sqrt((e[Sbest_interesting!=0.0]**2).mean())
        plt.loglog(my_time[fname], errors, color=my_color[fname], label=fname)
    else:
        print("We cannot compare with", fname, "because it doesn't have all the energies")
        print("  ", EminT,"<", my_energy[fname][0])
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'rms entropy error')

all_figures = set()
keep_going = True
while keep_going:
    keep_going = False
    for i in range(max_iter):
        for fig in all_figures:
            fig.clf()
        all_figures.add(plt.figure('Heat capacity'))
        plt.axvline(minT, linestyle=':', color='#ffaaaa')
        all_figures.add(plt.figure('Histogram'))
        plt.axvline(EminT, linestyle=':', color='#ffaaaa')

        all_figures.add(plt.figure('Normed entropy'))
        plt.axvline(EminT, linestyle=':', color='#ffaaaa')
        if args.match_energy is None:
            plt.plot(Ebest, Sbest - Sbest.max(), ':', color='#aaaaaa')
        else:
            plt.plot(Ebest, Sbest, ':', color='#aaaaaa')
        # all_figures.add(plt.figure('Temperature'))
        # plt.semilogy(Ebest_interesting,
        #              convex_hull_T(Ebest_interesting, Sbest_interesting), ':', color='#aaaaaa')
        for fname in args.yaml:
            if i < len(my_time[fname]):
                t = my_time[fname][i]
                j = i
            else:
                j = -1

            all_figures.add(plt.figure('Normed entropy'))
            if j > 0:
                if args.match_energy is None:
                    plt.plot(my_energy[fname],
                             my_entropy[fname][i-1,:]-my_entropy[fname][j-1,:].max(),
                             my_color[fname],
                             alpha=0.2)
                else:
                    plt.plot(my_energy[fname], my_entropy[fname][i-1,:], my_color[fname], alpha=0.2)
            if j == -1:
                if args.match_energy is None:
                    plt.plot(my_energy[fname],
                             my_entropy[fname][j,:]-my_entropy[fname][j,:].max(),
                             my_color[fname],
                             label=fname+' '+latex_float(len(my_entropy[fname])),
                             alpha=0.2)
                else:
                    plt.plot(my_energy[fname], my_entropy[fname][i-1,:], my_color[fname], alpha=0.2,
                             label=fname+' '+latex_float(len(my_entropy[fname])))
            else:
                if args.match_energy is None:
                    plt.plot(my_energy[fname],
                             my_entropy[fname][j,:]-my_entropy[fname][j,:].max(),
                             my_color[fname],
                             label=fname)
                else:
                    plt.plot(my_energy[fname], my_entropy[fname][i-1,:], my_color[fname], label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$S$')
            plt.legend(loc='best')
            # plt.ylim(Smin)

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

            all_figures.add(plt.figure('Heat capacity'))
            T = np.linspace(minT,1,1000)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('heat capacity')
            plt.xlabel('temperature')
            plt.plot(T, heat_capacity(T, my_energy[fname], my_entropy[fname][j,:]), my_color[fname],
                     label=fname)
            plt.legend(loc='best')
        plt.pause(0.000001)

plt.ioff()
plt.show()
