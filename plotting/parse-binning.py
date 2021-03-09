#!/usr/bin/python3

import yaml, sys, argparse, cbor, glob, itertools, os
import numpy as np
import scipy.constants as scipy
import readsystem

def latex_float(x):
    exp = int(np.log10(x*1.0))
    if abs(exp) > 2:
        x /= 10.0**exp
        if ('%.1g' % x) == '1':
            return r'10^{%.0f}' % (exp)
        return r'%.1g\times10^{%.0f}' % (x, exp)
    else:
        return '%g' % x

parser = argparse.ArgumentParser(description="parse data for energy histogram data")
parser.add_argument('yaml', nargs='*',
                    help = 'the names of some yaml files')
args = parser.parse_args()

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
            m = list(self.mean_extra('energy'))
            self._mean_energy = np.array([np.nan]+list(m)+[np.nan]) # pad with undefined energy in the unbounded bins
            assert(len(self._lnw) == len(self._hist))
        elif 'Linear' in data:
            self._kind = 'Linear'
            self._min = data['Linear']['min']
            self._width = data['Linear']['width']
            self._width = data['Linear']['width']
            self._lnw = np.array(data['Linear']['lnw']['total'])
            self._hist = np.array(data['Linear']['lnw']['count'])
            self._extra = data['Linear']['extra']
            m = list(self.mean_extra('energy'))
            self._mean_energy = np.array([np.nan]+list(m)+[np.nan]) # pad with undefined energy in the unbounded bins
            assert(len(self._lnw) == len(self._hist))
        else:
            self._kind = 'Histogram'
            print('bins has', data.keys())
            self._min = data['min']
            self._width = data['width']
            self._width = data['width']
            self._lnw = np.array(data['lnw'])
            self._hist = np.array(data['histogram'])
            self._extra = {}
            m = list(np.array(data['energy_total'])/np.array(data['histogram']))
            self._mean_energy = np.array([np.nan]+list(m)+[np.nan]) # pad with undefined energy in the unbounded bins
            assert(len(self._lnw) == len(self._hist))
        self._energy = np.arange(self._min + 0.5*self._width,
                                 self._min + len(self._lnw)*self._width,
                                 self._width)
        self._boundaries = np.arange(self._min,
                                     self._min + len(self._lnw)*self._width + 0.5*self._width,
                                     self._width)
    def excess_energy(self):
        if self._kind == 'Histogram':
            return self._energy
        elif self._kind == 'Linear':
            return self._energy
    def energy_boundaries(self):
        return self._boundaries
    def histogram(self, E=None):
        if self._kind == 'Histogram':
            return self._hist
        elif self._kind == 'Linear':
            return self._hist
    def lnw(self, E=None):
        if self._kind == 'Histogram':
            return self._lnw
        elif self._kind == 'Linear':
            return self._lnw
    def mean_extra(self, label):
        if self._kind == 'Histogram':
            return np.array(self._extra[label]['total'])/np.array(self._extra[label]['count'])
        elif self._kind == 'Linear':
            return np.array(self._extra[label]['total'])/np.array(self._extra[label]['count'])
    def extra_count(self, label):
        if self._kind == 'Histogram':
            return np.array(self._extra[label]['count'])
        elif self._kind == 'Linear':
            return np.array(self._extra[label]['count'])

class MC:
    """ The state of a Monte Carlo simulation """
    def __init__(self, filename):
        with open(filename, 'rb') as stream:
            if 'yaml' in filename:
                try:
                    self.data = yaml.full_load(stream)
                except IOError:
                    print('An error occurred trying to read the file.')
            elif 'cbor' in filename:
                try:
                    self.data = cbor.load(stream)
                except IOError:
                    print('An error occurred trying to read the file.')
            else:
                print('What kind of file is this?!')

            print('found', self.data.keys())
            self._bins = Bins(self.data['bins'])
            if 'high_resolution' in self.data and self.data['high_resolution'] is not None:
                self._high_resolution = Bins(self.data['high_resolution'])
            else:
                self._high_resolution = None
    def moves(self):
        return self.data['moves']
    def excess_energy(self):
        return self._bins.excess_energy()
    @property
    def energy_boundaries(self):
        return self._bins.energy_boundaries()
    @property
    def mean_energy(self):
        return self._bins._mean_energy
    @property
    def lnw(self):
        s = self.excess_entropy()
        s -= s.max()
        s -= np.log(np.sum(np.exp(s))) # normalize probability as best we can
        lnw = np.array([np.NINF]+list(s)+[np.NINF]) # pad with zero density of states in the unbounded bins
        # print('norm', np.sum(np.exp(lnw)))
        return lnw
    @property
    def system(self):
        return readsystem.readsystem(self.data)
    def histogram(self):
        return self._bins.histogram()
    def excess_entropy(self):
         lnw = 1.0*self._bins.lnw()
         if 'Sad' in self.data['method']:
             E = self.excess_energy()
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
             hist = self._bins.extra_count('hist')
             lnw[hist > 0] += np.log(hist[hist>0])
             if hist.sum() == 0:
                 print('where is the data?!')
                 print(hist)
         return lnw - lnw.max()

for fname in args.yaml:
    print(fname)
    mc = MC(fname)
    base = os.path.splitext(fname)[0]

    np.savetxt(base+'-energy-boundaries.dat', mc.energy_boundaries)
    np.savetxt(base+'-mean-energy.dat', mc.mean_energy)
    np.savetxt(base+'-lnw.dat', mc.lnw) #includes unbounded extremes

    with open(base+'-system.dat', 'w') as f:
        yaml.safe_dump(mc.system, f)