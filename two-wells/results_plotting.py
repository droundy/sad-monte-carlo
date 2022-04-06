import os
from typing import final
import numpy as np
import glob
import matplotlib.pyplot as plt
import system
import styles
import heat_capacity

class Results:
    def __init__(self):
        self.fnames = []
        self.mean_e=dict()
        self.mean_which=dict()
        self.hist=dict()
        self.E=dict()
        self.S=dict()
        self.T=dict()
        self.C=dict()
        self.moves=dict()
        self.errors_S=dict()
        self.errors_C=dict()
    
    def _query_fname(self, fname):
        data = dict()
        if fname not in self.fnames:
            return None
        if 'tem+' not in fname:
            data['mean_e']=self.mean_e[fname]
            data['mean_which']=self.mean_which[fname]
            data['hist']=self.hist[fname]
            data['E']=self.E[fname]
            data['S']=self.S[fname]
            data['errors_S']=self.errors_S[fname]

        data['T']=self.T[fname]
        data['C']=self.C[fname]
        data['moves']=self.moves[fname]
        data['errors_C']=self.errors_C[fname]
        return data
    
    def _mean_data(self, fname):
        min_data = np.nan_to_num(self._query_fname(fname), 
                                    nan=np.inf)
        max_data = np.nan_to_num(min_data.copy(), 
                                    nan=-np.inf)
        mean_data = np.nan_to_num(min_data.copy(),
                                    nan=0)

        given_seed = 'seed-'+fname.split('seed-')[-1].split('+')[0]
        #print(given_seed)
        filt = lambda f: \
            f.replace('seed-'+f.split('seed-')[-1].split('+')[0], '') == \
                fname.replace(given_seed, '') and f!=fname
        i=0
        for f in filter(filt, self.fnames):
            #print(f)
            i+=1
            new_data = self._query_fname(f)
            for k in new_data.keys():
                if new_data[k] is not None and k not in ['E','T']:
                    idx = min(len(min_data[k])-1, len(new_data[k])-1)
                    min_data[k] = np.minimum(min_data[k][:idx], np.nan_to_num(new_data[k][:idx], nan=np.inf))
                    max_data[k] = np.maximum(max_data[k][:idx], np.nan_to_num(new_data[k][:idx], nan=-np.inf))

                    mean_data[k] = mean_data[k][:idx] + np.nan_to_num(new_data[k][:idx], nan=0)
        for k in mean_data:
            if mean_data[k] is not None and k not in ['E','T']:
                mean_data[k] = mean_data[k] / i

        return min_data, mean_data, max_data


    def _stack_data_by_key(self, fname, key):
        data_stack = self._query_fname(fname)[key]

        given_seed = 'seed-'+fname.split('seed-')[-1].split('+')[0]
        #print(given_seed)
        filt = lambda f: \
            f.replace('seed-'+f.split('seed-')[-1].split('+')[0], '') == \
                fname.replace(given_seed, '') and f!=fname
        for f in filter(filt, self.fnames):
            #print(f)
            new_data = self._query_fname(f)[key]
            if new_data is not None:
                idx = min(data_stack.shape[-1], len(new_data))
                if len(data_stack.shape) == 1:
                    data_stack = np.vstack((data_stack[:idx], new_data[:idx]))
                else:
                    data_stack = np.vstack((data_stack[:,:idx], new_data[:idx]))
        return data_stack


    def _plot_from_data(self, ax, axins, fname, data=None, data_bounds=None, subplot = None, dump_into_thesis = None):
        if data is None:
            data = self._query_fname(fname)
        if data_bounds is not None:
            lower_data, upper_data = data_bounds
        base = fname[:-4]
        method = os.path.split(fname)[-1].split('+')[0]
        print(method)
        if method == 'itwl':
            label = r'$1/t$-WL' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[0]
        if method == 'sad':
            label = r'SAD' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[0]
        if method == 'z':
            label = r'ZMC' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[0]
        if method == 'tem':
            label = r'TEM' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[0]
        

        if method in {'wl','itwl','sad', 'z'}:
            plt.figure('fraction-well')
            plt.plot(data['mean_e'], data['mean_which'], label=label)
        
            if len(data['hist']) != 0:
                plt.figure('histogram')
                plt.plot(data['mean_e'], data['hist'], label=label)

            plt.figure('latest-entropy')
            plt.plot(data['E'][:len(data['S'])], data['S'], 
                                                    label=label, 
                                                    marker = styles.marker(base),
                                                    color = styles.color(base), 
                                                    linestyle= styles.linestyle(base), 
                                                    markevery=25)
            plt.figure('convergence')
            plt.loglog(data['moves'], data['errors_S'], 
                                                label=label, 
                                                marker = styles.marker(base), 
                                                color = styles.color(base), 
                                                linestyle= styles.linestyle(base), 
                                                markevery=4)
        elif method == 'z':
            plt.plot(data['E'], data['S'], 
                                    label=label, 
                                    color = styles.color(base), 
                                    linestyle= styles.linestyle(base))

            if data_bounds is not None:
                plt.fill_between(lower_data['moves'], lower_data['errors_S'], upper_data['errors_S'],
                                                    color = styles.color(base),
                                                    linestyle=styles.linestyle(base),
                                                    linewidth = 2,
                                                    alpha = 0.2)
        
        heat_capacity.plot_from_data(data['T'][:len(data['C'])], data['C'],
                                                                    fname=fname,
                                                                    ax=ax, 
                                                                    axins=axins)

        plt.figure('convergence-heat-capacity')
        plt.loglog(data['moves'], data['errors_C'], 
                                            label=label, 
                                            marker = styles.marker(base), 
                                            color = styles.color(base), 
                                            linestyle= styles.linestyle(base), 
                                            markevery=4)
        if data_bounds is not None:
            plt.fill_between(lower_data['moves'], lower_data['errors_C'], upper_data['errors_C'],
                                                color = styles.color(base),
                                                alpha = 0.2)
        
        if False:#subplot is not None:
            axs = subplot[0]
            axins_subplot = subplot[1]
            if method in {'wl','itwl','sad'}:
                axs['(c)'].plot(data['E'][:len(data['S'])], data['S'], 
                                                        label=label, 
                                                        marker = styles.marker(base),
                                                        color = styles.color(base), 
                                                        linestyle= styles.linestyle(base), 
                                                        markevery=250)
            elif method == 'z':
                axs['(c)'].plot(data['E'], data['S'], 
                                        label=label, 
                                        color = styles.color(base), 
                                        linestyle= styles.linestyle(base))
            
            heat_capacity.plot_from_data(data['T'][:len(data['C'])], data['C'],
                                                                        fname=fname,
                                                                        ax=axs['(d)'], 
                                                                        axins=axins_subplot)

            plt.figure('convergence')
            if method in {'wl','itwl','sad'}:
                axs['(a)'].loglog(data['moves'], data['errors_S'], 
                                                    label=label, 
                                                    marker = styles.marker(base), 
                                                    color = styles.color(base), 
                                                    linestyle= styles.linestyle(base), 
                                                    markevery=2)
            elif method == 'z':
                axs['(a)'].loglog(data['moves'], data['errors_S'], 
                                                    label=label, 
                                                    color = styles.color(base), 
                                                    linestyle= styles.linestyle(base), 
                                                    linewidth = 3)
            if data_bounds is not None:
                axs['(a)'].fill_between(lower_data['moves'], lower_data['errors_S'], upper_data['errors_S'],
                                                    color = styles.color(base),
                                                    linestyle=styles.linestyle(base),
                                                    linewidth = 2,
                                                    alpha = 0.2)


            plt.figure('convergence-heat-capacity')
            if method in {'wl','itwl','sad'}:
                axs['(b)'].loglog(data['moves'], data['errors_C'], 
                                                    label=label, 
                                                    marker = styles.marker(base), 
                                                    color = styles.color(base), 
                                                    linestyle= styles.linestyle(base), 
                                                    markevery=2)
            elif method == 'z':
                axs['(b)'].loglog(data['moves'], data['errors_C'], 
                                                    label=label, 
                                                    color = styles.color(base), 
                                                    linestyle= styles.linestyle(base), 
                                                    linewidth = 3)
            if data_bounds is not None:
                axs['(b)'].fill_between(lower_data['moves'], lower_data['errors_C'], upper_data['errors_C'],
                                                    color = styles.color(base),
                                                    alpha = 0.2)
        


    def add_npz(self, fname):
        if fname[-4:] != '.npz':
            raise('Incorrect filetype')
        self.fnames.append(fname)
        seed = fname.split('seed-')[-1].split('+')[0]
        data = np.load(fname)
        
        self.T[fname] = data['T']
        self.moves[fname] = data['moves']
        if 'tem+' not in fname:
            self.mean_e[fname] = data['mean_e']
            self.mean_which[fname] = data['mean_which']
            try:
                self.hist[fname] = data['hist']
            except:
                self.hist[fname] = None
            self.S[fname] = data['S']
            self.E[fname] = data['E']
            self.errors_S[fname] = data['errors_S']
        self.C[fname] = data['C']
        self.errors_C[fname] = data['errors_C']


    def plot_seed(self,
                    ax,
                    axins,
                    seed, 
                    method = None, 
                    additional_filters = None):
        seed = str(seed)
        filter_seed_hist = lambda f: 'seed-'+seed+'+' in f
        filter_seed_replicas = lambda f: 'seed-'+seed+'+' in f
        filters = [filter_seed_hist, filter_seed_replicas]
        if method is not None:
            method_filter = lambda f: method in f
        else:
            method_filter = lambda f: True
        filters.append(method_filter)
        if additional_filters is None:
            additional_filters = lambda f: True
        filters.append(additional_filters)
        fnames = self.fnames
        for filt in filters:
            fnames = filter(filt, fnames)

        for f in fnames:
            self._plot_from_data(ax, axins, f)
            

    def mean_method(self, ax, axins, subplot = None, dump_into_thesis = None):
        stacked_data = dict()
        unstacked_data = dict()
        for f in filter(lambda f: 'seed-1+' in f, self.fnames):
            for k in self._query_fname(f).keys():
                if 'error' in k or 'moves' in k: 
                    stacked_data[k] = self._stack_data_by_key(f, k)
                else:
                    unstacked_data[k] = self._query_fname(f)[k]
            min_data = dict()
            max_data = dict()
            for k in stacked_data.keys():
                min_data[k] = np.nanmin(stacked_data[k], axis=0)
                unstacked_data[k] = np.nanmean(stacked_data[k], axis=0)
                max_data[k] = np.nanmax(stacked_data[k], axis=0)
            #print(max_data['errors_S'])
            self._plot_from_data(ax, axins, f, 
                                    data = unstacked_data, 
                                    data_bounds=(min_data, max_data), 
                                    subplot = subplot, 
                                    dump_into_thesis = dump_into_thesis)
    
    def median_method(self, ax, axins, subplot = None, dump_into_thesis = None):
        stacked_data = dict()
        unstacked_data = dict()
        for f in filter(lambda f: 'seed-1+' in f, self.fnames):
            for k in self._query_fname(f).keys():
                if 'error' in k or 'moves' in k: 
                    stacked_data[k] = self._stack_data_by_key(f, k)
                else:
                    unstacked_data[k] = self._query_fname(f)[k]
            min_data = dict()
            max_data = dict()
            for k in stacked_data.keys():
                min_data[k] = np.nanmin(stacked_data[k], axis=0)
                unstacked_data[k] = np.nanmedian(stacked_data[k], axis=0)
                max_data[k] = np.nanmax(stacked_data[k], axis=0)
            #print(max_data['errors_S'])
            self._plot_from_data(ax, axins, f, 
                                    data = unstacked_data, 
                                    data_bounds=(min_data, max_data), 
                                    subplot = subplot, 
                                    dump_into_thesis = dump_into_thesis)
