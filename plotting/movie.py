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
fnames = sys.argv[1:]
for fname in fnames:
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
print('energies:', Ebest[ind_minT], Ebest[ind_maxS])
Sbest_interesting = Sbest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]
Ebest_interesting = Ebest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]

plt.ion()

plt.figure('gamma')
for fname in fnames:
        plt.loglog(my_gamma_t[fname], my_gamma[fname], color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'$\gamma$')
plt.ylim(1e-12, 1.1)

def convex_hull(S):
    convexS = np.zeros_like(S)
    if len(convexS) > 1000:
        convexS[:] = S
        return convexS
    for i in range(len(S)):
        if S[i] > 0 and S[i] >= convexS[i]:
            for j in range(i+1,len(S)):
                if S[j] > 0 and S[j] >= convexS[j]:
                    for k in range(i, j+1):
                        convexS[k] = max(convexS[k], (S[i]*(j-k) + S[j]*(k-i))/(j-i))
    return convexS
def convex_hull_T(E, S):
    convexS = convex_hull(S)
    T = np.zeros_like(S)
    if convexS[1] > convexS[0]:
        T[0] = (E[1]-E[0])/(convexS[1]-convexS[0])
    if convexS[-1] > convexS[-2]:
        T[-1] = (E[-1]-E[-2])/(convexS[-1]-convexS[-2])
    for i in range(1,len(T)-1):
        if convexS[i+1] > convexS[i-1]:
            T[i] = (E[i+1]-E[i-1])/(convexS[i+1]-convexS[i-1])
    return T

Tbest_interesting = convex_hull_T(Ebest_interesting, Sbest_interesting)
plt.figure('temperature-comparison')
for fname in my_energy.keys():
    if my_energy[fname][0] <= EminT:
        errors = np.zeros(len(my_time[fname]))
        ind_minT = np.argwhere(my_energy[fname] == EminT)[0][0]
        ind_maxS = np.argwhere(my_energy[fname] == EmaxS)[0][0]
        E_interesting = my_energy[fname][ind_minT:ind_maxS+1]
        for i in range(len(my_time[fname])):
            S_interesting = my_entropy[fname][i,ind_minT:ind_maxS+1]
            T_interesting = convex_hull_T(E_interesting, S_interesting)
            e = (T_interesting - Tbest_interesting)/Tbest_interesting
            errors[i] = np.sqrt((e**2).mean())
        plt.loglog(my_time[fname], errors, color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'rms relative error in $T$')

plt.figure('comparison')
for fname in fnames:
    if my_energy[fname][0] <= EminT:
        errors = np.zeros(len(my_time[fname]))
        ind_minT = np.argwhere(my_energy[fname] == EminT)[0][0]
        ind_maxS = np.argwhere(my_energy[fname] == EmaxS)[0][0]
        for i in range(len(my_time[fname])):
            S_interesting = my_entropy[fname][i,ind_minT:ind_maxS+1]
            e = S_interesting - Sbest_interesting
            e -= e.mean()
            errors[i] = np.sqrt((e**2).mean())
        plt.loglog(my_time[fname], errors, color=my_color[fname], label=fname)
    else:
        print("We cannot compare with", fname, "because it doesn't have all the energies")
        print("  ", EminT,"<", my_energy[fname][0])
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
        all_figures.add(plt.figure('Temperature'))
        plt.semilogy(Ebest_interesting,
                     convex_hull_T(Ebest_interesting, Sbest_interesting), ':', color='#aaaaaa')
        for fname in fnames:
            if i < len(my_time[fname]):
                t = my_time[fname][i]
                j = i
            else:
                j = -1
            if fname == fnames[0]:
                print('frame', i, 'with', t, 'iterations')

            # all_figures.add(plt.figure('Entropy'))
            # if j > 0:
            #     plt.plot(my_energy[fname], my_entropy[fname][j-1,:], my_color[fname],
            #              alpha=0.2)
            # if j == -1:
            #     plt.plot(my_energy[fname], my_entropy[fname][j,:], my_color[fname],
            #              label=fname+' '+latex_float(len(my_entropy[fname])),
            #              alpha=0.2)
            # else:
            #     plt.plot(my_energy[fname], my_entropy[fname][j,:], my_color[fname],
            #              label=fname)
            # plt.title('$t=%s/%s$' % (latex_float(t),
            #                          latex_float(my_time[fname][-1])))
            # plt.ylabel('$S$')
            # plt.legend(loc='best')

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
                # plt.plot(my_energy[fname],
                #          convex_hull(my_entropy[fname][j,:])-my_entropy[fname][j,:].max(),
                #          ':',
                #          color=my_color[fname],
                #          label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$S$')
            plt.legend(loc='best')
            plt.ylim(Smin, 0)

            all_figures.add(plt.figure('Temperature'))
            T = convex_hull_T(my_energy[fname], my_entropy[fname][j,:])
            if len(T[T>0]) > 1:
                if j == -1:
                    plt.semilogy(my_energy[fname][T>0],
                                 T[T>0],
                                 color=my_color[fname],
                                 label=fname+' '+latex_float(len(my_entropy[fname])),
                                 alpha=0.5)
                else:
                    plt.semilogy(my_energy[fname][T>0],
                                 T[T>0],
                                 color=my_color[fname],
                                 label=fname)
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('$T$')
            plt.legend(loc='best')
            # plt.ylim(Tbest_interesting.min(), Tbest_interesting.max())

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
