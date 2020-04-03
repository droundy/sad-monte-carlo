#!/usr/bin/python3

import yaml, sys, argparse, cbor, glob
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

class MC:
    """ The state of a Monte Carlo simulation """
    def __init__(self, filename):
        with open(filename, 'rb') as f:
            self.data = cbor.load(f)
            self._bins = Bins(self.data['bins'])
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
         lnw = self._bins.lnw()
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
         return lnw

plt.ion()
assert(len(args.dirname)==1) # fix this later
for f in sorted(glob.glob(args.dirname[0]+'/*.cbor')):
    mc = MC(f)

    plt.figure('histogram')
    plt.clf()
    plt.title(f)
    plt.plot(mc.energy(), mc.histogram())

    plt.figure('lnw')
    plt.clf()
    plt.title(f)
    plt.plot(mc.energy(), mc.entropy())

    # plt.show()
    plt.pause(0.1)
plt.ioff()
plt.show()
