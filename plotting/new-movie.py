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
Smin = None
minT = 0.01

def fix_fname(fname):
    if fname[-5:] == '.yaml':
        return fname[:-5]
    return fname

fnames = [fix_fname(f) for f in sys.argv[1:]]
for fname in fnames:
    print(fname)
    my_histogram[fname] = np.loadtxt(fname+'.histogram')
    my_energy[fname] = np.loadtxt(fname+'.energy')
    my_entropy[fname] = np.loadtxt(fname+'.entropy')
    my_time[fname] = np.loadtxt(fname+'.time')
    my_color[fname] = allcolors.pop()
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    if Smin is None:
        Ebest = my_energy[fname];
        Sbest = my_entropy[fname][-1,:]
        Smin = Sbest[Sbest!=0].min() - Sbest.max()

EmaxS = Ebest[np.argmax(Sbest)]
EminT = Ebest[np.argmax(Sbest*minT - Ebest)]
ind_minT = np.argwhere(Ebest == EminT)[0][0]
ind_maxS = np.argwhere(Ebest == EmaxS)[0][0]
print('energies:', Ebest[ind_minT], 'at temperature', minT, 'and max entropy', Ebest[ind_maxS])
Sbest_interesting = Sbest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]
Ebest_interesting = Ebest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]

plt.ion()

def heat_capacity(T, E, S):
    C = np.zeros_like(T)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P = P/P.sum()
        U = (E*P).sum()
        C[i] = ((E-U)**2*P).sum()/T[i]**2
    return C

# Tbest_interesting = convex_hull_T(Ebest_interesting, Sbest_interesting)
# plt.figure('temperature-comparison')
# for fname in my_energy.keys():
#     if my_energy[fname][0] <= EminT:
#         errors = np.zeros(len(my_time[fname]))
#         ind_minT = np.argwhere(my_energy[fname] == EminT)[0][0]
#         ind_maxS = np.argwhere(my_energy[fname] == EmaxS)[0][0]
#         E_interesting = my_energy[fname][ind_minT:ind_maxS+1]
#         for i in range(len(my_time[fname])):
#             S_interesting = my_entropy[fname][i,ind_minT:ind_maxS+1]
#             T_interesting = convex_hull_T(E_interesting, S_interesting)
#             e = (T_interesting - Tbest_interesting)/Tbest_interesting
#             errors[i] = np.sqrt((e**2).mean())
#         plt.loglog(my_time[fname], errors, color=my_color[fname], label=fname)
# plt.legend(loc='best')
# plt.xlabel('$t$')
# plt.ylabel(r'rms relative error in $T$')

def indwhere(array, value):
    best = 0
    for i in range(len(array)):
        if abs(array[i]-value) < abs(array[best]-value):
            best = i
    return best

def Sbest_function(e):
    return np.interp(e, Ebest_interesting, Sbest_interesting)

plt.figure('comparison')
for fname in fnames:
    if my_energy[fname][0] <= EminT:
        errors = np.zeros(len(my_time[fname]))
        ind_minT = indwhere(my_energy[fname], EminT)
        ind_maxS = indwhere(my_energy[fname], EmaxS)
        for i in range(len(my_time[fname])):
            S_interesting = my_entropy[fname][i,ind_minT:ind_maxS+1]
            E_interesting = my_energy[fname][ind_minT:ind_maxS+1]
            e = np.zeros_like(S_interesting)
            for j in range(len(e)):
                e[j] = S_interesting[j] - Sbest_function(E_interesting[j])
            e -= e.mean() # WARNING, this causes problems if there are impossible states in the interesting energy range.
            errors[i] = np.sqrt((e**2).mean()) # WARNING, this causes problems if there are impossible states in the interesting energy range.
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
        for E0 in np.linspace(2*Ebest.min() - Ebest.max(), Ebest.max(), 20):
            plt.plot(Ebest, Smin + (Ebest - E0)/minT, ':', color='#ffeedd')
        plt.axvline(EminT, linestyle=':', color='#ffaaaa')
        plt.plot(Ebest, Sbest - Sbest.max(), ':', color='#aaaaaa')
        # all_figures.add(plt.figure('Temperature'))
        # plt.semilogy(Ebest_interesting,
        #              convex_hull_T(Ebest_interesting, Sbest_interesting), ':', color='#aaaaaa')
        for fname in fnames:
            if i < len(my_time[fname]):
                t = my_time[fname][i]
                j = i
            else:
                j = -1
            #if fname == fnames[0]:
            #    print('frame', i, 'with', t, 'iterations')

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

            # all_figures.add(plt.figure('Temperature'))
            # T = convex_hull_T(my_energy[fname], my_entropy[fname][j,:])
            # if len(T[T>0]) > 1:
            #     if j == -1:
            #         plt.semilogy(my_energy[fname][T>0],
            #                      T[T>0],
            #                      color=my_color[fname],
            #                      label=fname+' '+latex_float(len(my_entropy[fname])),
            #                      alpha=0.5)
            #     else:
            #         plt.semilogy(my_energy[fname][T>0],
            #                      T[T>0],
            #                      color=my_color[fname],
            #                      label=fname)
            # plt.title('$t=%s/%s$' % (latex_float(t),
            #                          latex_float(my_time[fname][-1])))
            # plt.ylabel('$T$')
            # plt.legend(loc='best')
            # # plt.ylim(Tbest_interesting.min(), Tbest_interesting.max())

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
        plt.pause(0.1)

plt.ioff()
plt.show()
