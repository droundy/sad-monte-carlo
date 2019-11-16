#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import re, argparse, os
import martiniani

parser = argparse.ArgumentParser(description="create movie and graphs for lj-cluster data")
parser.add_argument('--minT', action='store', type=float, default=0.005,
                    help = "the minimum temperature of interest")
parser.add_argument('--match-energy', action='store', type=float,
                    help = "the energy at which we want to normalize the entropy")
parser.add_argument('yaml', nargs='*',
                    help = 'the names of some yaml files')
args = parser.parse_args()

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


data4 = np.loadtxt('LJ31_Cv_Reference_4.csv', delimiter = ',', unpack = True)
T = data4[0]
CV = (data4[1]-3/2)*31

T = np.linspace(args.minT, 0.4, 1000)

other_T = martiniani.T
other_CV = martiniani.CV
other_name = 'Martiniani et al.'
j_lower_peak = 0
for j in range(len(T)):
    if T[j] < 0.2:
        j_lower_peak = j

my_energy = {}
my_de = {}
my_histogram = {}
my_entropy = {}
my_time = {}
my_color = {}
max_iter = 0
Smin = None
Smax = 0
my_minT = {}
my_too_lo = {}
my_too_hi = {}
minT = args.minT

my_cv_error = {}
my_cv_peak_T = {}
my_cv_peak = {}
cv_iters = []

def lookup_entry(entry, yaml_data):
    x = re.search('{}: (\S+)'.format(entry), yaml_data)
    if x:
        return float(x.group(1))

def fix_fname(fname):
    if fname[-5:] == '.yaml':
        return fname[:-5]
    if fname[-5:] == '.cbor':
        return fname[:-5]
    return fname

print(args)
fnames = [fix_fname(f) for f in args.yaml]

def heat_capacity(T, E, S):
    C = np.zeros_like(T)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P = P/P.sum()
        U = (E*P).sum()
        C[i] = ((E-U)**2*P).sum()/T[i]**2
    return C

for fname in fnames:
    print(fname)
    if os.path.exists(fname+'.yaml'):
        with open(fname+'.yaml') as f:
            yaml = f.read()
            my_too_hi[fname] = lookup_entry('too_hi', yaml)
            my_too_lo[fname] = lookup_entry('too_lo', yaml)
            my_minT[fname] = lookup_entry('min_T', yaml)
    else:
        my_too_lo[fname] = None
        my_too_hi[fname] = None
        my_minT[fname] = None
    my_time[fname] = np.loadtxt(fname+'.time')
    first_frame = 0
    for i in range(len(my_time[fname])):
        if my_time[fname][i] <= 1e8:
            first_frame = i
    my_time[fname] = my_time[fname][first_frame:]
    my_histogram[fname] = np.loadtxt(fname+'.histogram')[first_frame:]
    my_energy[fname] = np.loadtxt(fname+'.energy')
    my_de[fname] = my_energy[fname][1] - my_energy[fname][0]
    my_entropy[fname] = np.loadtxt(fname+'.entropy')[first_frame:]
    my_cv_error[fname] = []
    my_cv_peak[fname] = []
    my_cv_peak_T[fname] = []
    if minT is None and my_minT[fname] is not None:
        minT = my_minT[fname]
    my_color[fname] = allcolors.pop()
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    norm_entropy = lambda s: s.max()
    if args.match_energy is not None:
        match_index = (np.abs(my_energy[fname] - args.match_energy)).argmin()
        print('match_index is', match_index, 'from', args.match_energy)
        norm_entropy = lambda s: s[match_index]
    if Smin is None:
        Ebest = my_energy[fname];
        Sbest = my_entropy[fname][-1,:] - norm_entropy(my_entropy[fname][-1,:])
        CV = heat_capacity(T, Ebest, Sbest)
        np.savetxt("best_cv.txt", np.array([T, CV]).transpose())
        Smin = Sbest[Sbest!=0].min()
    Smax = max(Smax, (my_entropy[fname][-1,:] - norm_entropy(my_entropy[fname][-1,:])).max())
    print('Smax is now', Smax)

EmaxS = Ebest[np.argmax(Sbest)]
EminT = Ebest[np.argmax(Sbest*minT - Ebest)]
ind_minT = np.argwhere(Ebest == EminT)[0][0]
ind_maxS = np.argwhere(Ebest == EmaxS)[0][0]
print('energies:', Ebest[ind_minT], 'at temperature', minT, 'and max entropy', Ebest[ind_maxS])
Sbest_interesting = Sbest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]
Ebest_interesting = Ebest[np.argwhere(Ebest == EminT)[0][0]:np.argwhere(Ebest == EmaxS)[0][0]+1]

plt.ion()

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

            norm_entropy = lambda s: s.max()
            if args.match_energy is not None:
                match_index = (np.abs(my_energy[fname] - args.match_energy)).argmin()
                norm_entropy = lambda s: s[match_index]
            all_figures.add(plt.figure('Normed entropy'))
            if my_too_lo[fname]:
                plt.axvline(my_too_lo[fname], linestyle=':', color=my_color[fname])
            if j > 0:
                plt.plot(my_energy[fname],
                         my_entropy[fname][i-1,:]-norm_entropy(my_entropy[fname][j-1,:]),
                         my_color[fname],
                         alpha=0.2)
            if j == -1:
                plt.plot(my_energy[fname],
                         my_entropy[fname][j,:]-norm_entropy(my_entropy[fname][j,:]),
                         my_color[fname],
                         label=fname+' '+latex_float(len(my_entropy[fname])),
                         alpha=0.2)
            elif my_entropy[fname].shape[0] > j and my_entropy[fname].shape[1] > j:
                plt.plot(my_energy[fname],
                         my_entropy[fname][j,:]- norm_entropy(my_entropy[fname][j,:]),
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
            plt.ylim(Smin, Smax)

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
            if my_too_lo[fname]:
                plt.axvline(my_too_lo[fname], linestyle=':', color=my_color[fname])
            if my_too_hi[fname]:
                plt.axvline(my_too_hi[fname], linestyle=':', color=my_color[fname])
            plt.ylabel('histogram')
            if j > 0:
                plt.plot(my_energy[fname], my_histogram[fname][j-1,:]/my_de[fname], my_color[fname],
                         alpha=0.2)
            if j == -1:
                plt.plot(my_energy[fname], my_histogram[fname][j,:]/my_de[fname], my_color[fname],
                         label=fname+' '+latex_float(len(my_entropy[fname])),
                         alpha=0.2)
            else:
                plt.plot(my_energy[fname], my_histogram[fname][j,:]/my_de[fname], my_color[fname],
                         label=fname)
            plt.legend(loc='best')

            all_figures.add(plt.figure('Heat capacity'))
            plt.title('$t=%s/%s$' % (latex_float(t),
                                     latex_float(my_time[fname][-1])))
            plt.ylabel('heat capacity')
            plt.xlabel('temperature')
            mycv = heat_capacity(T, my_energy[fname], my_entropy[fname][j,:])
            plt.plot(T, mycv, my_color[fname], label=fname)
            plt.legend(loc='best')

            err = 0
            norm = 0
            peak = 0
            peak_T = 0
            for j in range(1,j_lower_peak):
                err += (T[j+1]-T[j-1])*abs(CV[j]-mycv[j])
                norm += (T[j+1]-T[j-1])
                if mycv[j] > peak:
                    peak = mycv[j]
                    peak_T = T[j]
            if len(my_cv_error[fname]) < len(my_time[fname]):
                my_cv_error[fname].append(err/norm)
                my_cv_peak[fname].append(peak)
                my_cv_peak_T[fname].append(peak_T)

            all_figures.add(plt.figure('CV errors'))
            plt.loglog(my_time[fname][:len(my_cv_error[fname])], my_cv_error[fname],
                           my_color[fname], label=fname)
            plt.xlabel('number of moves')
            plt.ylabel('mean error in $C_V$')
            plt.legend(loc='best')

            true_cv_peak = 75.645
            true_cv_peak_T = 0.02745

            # all_figures.add(plt.figure('CV peak temperature'))
            # plt.semilogx(my_time[fname][:len(my_cv_peak_T[fname])],
            #              np.array(my_cv_peak_T[fname]),
            #              my_color[fname], label=fname)
            # plt.xlabel('number of moves')
            # plt.ylabel('error in peak temp')
            # plt.axhline(true_cv_peak_T, color='k', linewidth=0.01)
            # plt.ylim(0,4*true_cv_peak_T)
            # plt.legend(loc='best')

            # all_figures.add(plt.figure('CV peak value'))
            # plt.semilogx(my_time[fname][:len(my_cv_peak[fname])],
            #              np.array(my_cv_peak[fname]),
            #              my_color[fname], label=fname)
            # plt.xlabel('number of moves')
            # plt.ylabel('error in peak value')
            # plt.axhline(true_cv_peak, color='k', linewidth=0.01)
            # plt.ylim(0,10*true_cv_peak)
            # plt.legend(loc='best')

        all_figures.add(plt.figure('Heat capacity'))
        plt.plot(T, CV, 'k:', label='ref 4?')
        plt.plot(other_T, other_CV, 'k-', label=other_name)
        plt.ylim(0,140)
        plt.legend(loc='best')

        plt.pause(0.1)

plt.ioff()
plt.show()
