#!/usr/bin/python3

import yaml, sys, argparse, cbor, glob, itertools
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
parser.add_argument('dirname', nargs='*',
                    help = 'the names of some yaml files')
args = parser.parse_args()

allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
                           'xkcd:lightblue', 'xkcd:puke', 'xkcd:puce', 'xkcd:turquoise']))

class Bins:
    """ A thing with bins """
    def __init__(self, data):
        if 'Histogram' in data:
            self._kind = 'Histogram'
            self._min = data['Histogram']['min']
            self._width = data['Histogram']['width']
            self._width = data['Histogram']['width']
            self._lnw = np.array(data['Histogram']['lnw']['total'])
            self._hist = np.array(data['Histogram']['lnw']['count'])
            self._extra = data['Histogram']['extra']
            assert(len(self._lnw) == len(self._hist))
            self._energy = np.arange(self._min + 0.5*self._width,
                                     self._min + len(self._lnw)*self._width,
                                     self._width)
    def energy(self):
        if self._kind == 'Histogram':
            return self._energy
    def histogram(self, E=None):
        if self._kind == 'Histogram':
            return self._hist
    def lnw(self, E=None):
        if self._kind == 'Histogram':
            return self._lnw
    def mean_extra(self, label):
        if self._kind == 'Histogram':
            return np.array(self._extra[label]['total'])/np.array(self._extra[label]['count'])

class MC:
    """ The state of a Monte Carlo simulation """
    def __init__(self, filename):
        with open(filename, 'rb') as f:
            self.data = cbor.load(f)
            self._bins = Bins(self.data['bins'])
    def moves(self):
        return self.data['moves']
    def energy(self):
        return self.data['system']['E']
    def cell(self):
        diag = self.data['system']['cell']['box_diagonal']
        return np.array([diag['x'], diag['y'], diag['z']])
    def volume(self):
        cell = self.cell()
        return cell[0]*cell[1]*cell[2]
    def energy(self):
        return self._bins.energy()
    def histogram(self):
        return self._bins.histogram()
    def entropy(self):
         lnw = 1.0*self._bins.lnw()
         if 'Sad' in self.data['method']:
             E = self.energy()
             hist = self.histogram()
             too_lo = self.data['method']['Sad']['too_lo']
             too_hi = self.data['method']['Sad']['too_hi']
             i_lo = abs(E - too_lo).argmin()
             i_hi = abs(E - too_hi).argmin()
             min_T = self.data['method']['Sad']['min_T']
             mean_hist = hist[i_lo:i_hi+1].mean()
             lnw[E < too_lo] = lnw[i_lo] + (E[E<too_lo] - too_lo)/min_T + np.log(hist[E<too_lo]/mean_hist)
             lnw[E > too_hi] = lnw[i_hi] + np.log(hist[E>too_hi]/mean_hist)
             lnw[np.isnan(lnw)] = 0
             lnw[np.isinf(lnw)] = 0
         elif 'WL' in self.data['method'] and self.data['method']['WL']['gamma'] == 0:
             hist = np.array(self.data['bins']['Histogram']['extra']['hist']['count'])
             lnw[hist > 0] += np.log(hist[hist>0])
             if hist.sum() == 0:
                 print('where is the data?!')
                 print(hist)
         return lnw
    def temperature(self):
        energy = self.energy()
        entropy = self.entropy()
        temp = np.zeros_like(entropy)
        if len(entropy) < 2:
            return temp
        for i in range(0, len(energy)):
            if i == 0:
                dE = energy[i+1] - energy[i]
                dS = entropy[i+1] - entropy[i]
            elif i == len(energy)-1:
                dE = energy[i] - energy[i-1]
                dS = entropy[i] - entropy[i-1]
            else:
                dE = energy[i+1] - energy[i-1]
                dS = entropy[i+1] - entropy[i-1]
            try:
                temp[i] = float(dE/dS)
            except ZeroDivisionError:
                if dE > 0:
                    temp[i] = -float('inf')
                temp[i] = float('inf') #if dE <= 0
        return temp
    def excess_pressure(self):
        if 'pressure' not in self._bins._extra:
            return np.zeros_like(self.energy())
        return self._bins.mean_extra('pressure')
    # TO DO: add temperature method, add pressure method (any more?)

    # We can also compute chemical potential mu_exc (and later mu) from:
    # G = U - TS + pV
    # G_exc = U_exc - T S_exc + p_exc V
    # mu_exc = G_exc/N
    def find_entropy(self, E):
         lnw = self.entropy()
         e = self.energy()
         return np.interp(E, e, lnw, left=0, right=0)

plt.ion()
assert(len(args.dirname)>=1) # we need at least one dirname
reference = sorted(glob.glob(args.dirname[0]+'/*.cbor'))[-1]
ref = MC(reference)

def max_entropy_error(mc):
    E = ref.energy()
    good_S = ref.entropy()
    bad_S = mc.find_entropy(E)
    return (good_S - bad_S).max() - (good_S - bad_S).min()

things = [sorted(glob.glob(dirname+'/*.cbor')) for dirname in args.dirname]
things = itertools.zip_longest(*things)
all_moves = []
all_max_entropy_errors = {}
all_last = {}
all_labels = {}
all_figures = set({})
for fs in things:
    for fig in all_figures:
        fig.clf()
    plt.figure('entropy')
    plt.plot(ref.energy(), ref.entropy() - ref.entropy().max(), ':', color='gray', label=reference)
    added_move = False
    for i in range(len(fs)):

        if fs[i] is None:
            mc = all_last[i]
            alpha = 0.5
            label = all_labels[i]
        else:
            f = fs[i]
            mc = MC(f)
            alpha = 1.0

            label = f.split('/')[0]
            all_labels[i] = label
            moves = mc.moves() # float(int(f.split('/')[1].split('.')[0]))
            title = '${}$'.format(latex_float(moves))

        all_last[i] = mc

        if not added_move:
            added_move = True
            all_moves.append(moves)
        if label not in all_max_entropy_errors:
            all_max_entropy_errors[label] = []
        all_max_entropy_errors[label].append(max_entropy_error(mc))

        all_figures.add(plt.figure('histogram'))
        plt.title(title)
        plt.plot(mc.energy(), mc.histogram(), label=label, alpha=alpha)
        plt.xlabel('$E$')
        plt.ylabel('$H$')
        plt.legend(loc='best')

        all_figures.add(plt.figure('entropy'))
        plt.title(title)
        plt.plot(mc.energy(), mc.entropy() - mc.entropy().max(), label=label, alpha=alpha)
        plt.xlabel('$E$')
        plt.ylabel('$S_{exc}$')
        plt.legend(loc='best')

        all_figures.add(plt.figure('energy_temperature'))
        plt.title(title)
        plt.plot(mc.energy(), mc.temperature(), label=label, alpha=alpha)
        plt.xlabel('$E$')
        plt.ylabel('$T$')
        plt.legend(loc='upper right')

        all_figures.add(plt.figure('excess pressure'))
        plt.title(title)
        plt.plot(mc.energy(), mc.excess_pressure(), label=label, alpha=alpha)
        plt.xlabel('$E$')
        plt.ylabel('$p_{exc}$')
        plt.legend(loc='best')

    # plt.show()
    plt.pause(0.1)

plt.figure('comparison')
for label in all_labels.values():
    plt.loglog(all_moves, all_max_entropy_errors[label], '.-', label=label)
plt.legend(loc='best')

plt.ioff()
plt.show()
