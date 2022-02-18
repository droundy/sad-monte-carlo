#!/usr/bin/python3
import os
import numpy as np
import matplotlib.pyplot as plt
import system, styles
import find_phase_transition
import time
from mytimer import Timer

T_peak = find_phase_transition.actual_T

def C(T, S):#T is a temperature and S is an entropy function
    # start = time.process_time()
    E = np.linspace(-system.h_small, 0, 1000)
    E = 0.5*(E[1:] + E[:-1])
    dE = E[1] - E[0]

    def normalize_S(S):
        S = S - max(S)
        total = np.sum(np.exp(S)*dE)
        return S - np.log(total)

    S = S(E) # normalize_S(S(E))

    S_minus_E = S-E/T
    M = np.max(S_minus_E)

    Z = np.sum(np.exp(S_minus_E-M))*dE

    avg_E = np.sum(np.exp(S_minus_E-M) * E) * dE / Z

    avg_E_squared = np.sum(np.exp(S_minus_E-M) * E**2) * dE / Z

    # print('C took', time.process_time() - start)
    return (avg_E_squared - avg_E**2 ) / T**2


def C_E_Esqrd(T, S):#T is a temperature and S is an entropy function
    # start = time.process_time()
    E = np.linspace(-system.h_small, 0, 1000)
    E = 0.5*(E[1:] + E[:-1])
    dE = E[1] - E[0]

    def normalize_S(S):
        S = S - max(S)
        total = np.sum(np.exp(S)*dE)
        return S - np.log(total)

    S = S(E) # normalize_S(S(E))

    S_minus_E = S-E/T
    M = np.max(S_minus_E)

    Z = np.sum(np.exp(S_minus_E-M))*dE

    avg_E = np.sum(np.exp(S_minus_E-M) * E) * dE / Z

    avg_E_squared = np.sum(np.exp(S_minus_E-M) * E**2) * dE / Z

    # print('C took', time.process_time() - start)
    return ((avg_E_squared - avg_E**2 ) / T**2, avg_E, avg_E_squared)
def C_vector(T, S):#T is an array of temperatures and S is an entropy function
    # start = time.process_time()
    E = np.linspace(-system.h_small, 0, 1000)
    E = 0.5*(E[1:] + E[:-1])
    dE = E[1] - E[0]

    def normalize_S(S):
        S = S - max(S)
        total = np.sum(np.exp(S)*dE)
        return S - np.log(total)

    S = S(E) # normalize_S(S(E))
    E_integrate = np.tile(E, (len(T),1))

    S_minus_E = np.array([S-E/t for t in T])
    M = np.max(S_minus_E, axis=1)
    print(M.shape)

    Z = np.sum(np.exp(S_minus_E-M), axis=1)*dE

    avg_E = np.sum(np.exp(S_minus_E-M) * E, axis=1) * dE / Z

    avg_E_squared = np.sum(np.exp(S_minus_E-M) * E**2, axis=1) * dE / Z

    # print('C took', time.process_time() - start)
    return (avg_E_squared - avg_E**2 ) / T**2


def _set_temperatures(ax=None, axins=None, Tmax=0.25):
    T_width = T_peak/2 # this is just a guess
    t_low = np.linspace(T_peak/10,T_peak - T_width,10)
    t_peak = np.linspace(T_peak - T_width,T_peak + T_width,150)
    if axins is not None:
        axins.set_xlim(0, T_peak + T_width+0.005)
        axins.set_ylim(5, 27)
    t_high = np.linspace(T_peak + T_width,Tmax, 200)
    if ax is not None:
        ax.set_xlim(0,max(t_high))
        ax.set_ylim(0, 30)
    return (t_low, t_peak, t_high)


def plot(S, fname=None, ax=None, axins=None, Tmax=0.25):
    timer = Timer(f'heat_capacity.plot {fname}')
    if fname is not None:
        base = fname[:-8]
        method = base[:base.find('-')]
    else:
        base = None
        method = None

    t_low, t_peak, t_high = _set_temperatures(ax=ax,axins=axins,Tmax=Tmax)

    timer_low_T = Timer('C for low T')
    c_low = np.array([C(T,S) for T in t_low])
    del timer_low_T

    timer_peak = Timer('C for peak T')
    c_peak = np.array([C(T,S) for T in t_peak])
    del timer_peak

    c_high = np.array([C(T,S) for T in t_high])

    ax.plot(np.concatenate((t_low,t_peak,t_high)), 
            np.concatenate((c_low,c_peak,c_high)), 
            label=styles.pretty_label(base), 
            color = styles.color(base), 
            linestyle= styles.linestyle(base), 
            markevery=10)

    # inset axes....
    #axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])
    axins.plot(np.concatenate((t_low,t_peak,t_high)), 
               np.concatenate((c_low,c_peak,c_high)), 
               label=styles.pretty_label(base), 
               marker = styles.marker(base), 
               color = styles.color(base), 
               linestyle= styles.linestyle(base), 
               markevery=10)
    # sub region of the original image
    # x1, x2, y1, y2 = 0.002, 0.009, 25, 140
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_xticklabels('')
    # axins.set_yticklabels('')

    # ax.indicate_inset_zoom(axins, edgecolor="black")




def plot_from_data(T_data, C_data, fname=None, ax=None, axins=None, Tmax=0.25):
    if fname is not None:
        base = fname[:-4]
        method = base[:base.find('-')]
    else:
        base = None
        method = None

    _, t_peak, _ = _set_temperatures(ax=ax,axins=axins,Tmax=Tmax)

    ax.plot(T_data, 
            C_data, 
            label=base,
            marker = styles.marker(base),
            color = styles.color(base), 
            linestyle= styles.linestyle(base), 
            markevery=10)

    # inset axes....
    #axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])
    axins.plot(T_data, 
               C_data, 
               label=styles.pretty_label(base), 
               marker = styles.marker(base), 
               color = styles.color(base), 
               linestyle= styles.linestyle(base), 
               markevery=10)
    x1, x2, y1, y2 = 0.002, 0.009, 5, 30
    # axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)



# returns data for heat capcity plot. Only returns single 
# temperature and heat capacity arrays
def data(S, fname=None, Tmax=0.25):
    # T_width = T_peak/2 # this is just a guess
    # t_low = np.linspace(T_peak/10,T_peak - T_width,10)
    # t_peak = np.linspace(T_peak - T_width,T_peak + T_width,150)
    # t_high = np.linspace(T_peak + T_width,Tmax, 10)
    t_low, t_peak, t_high = _set_temperatures()
    t = np.concatenate( (t_low,t_peak,t_high) )
    

    c_low = np.array([C(T,S) for T in t_low])
    c_peak = np.array([C(T,S) for T in t_peak])
    c_high = np.array([C(T,S) for T in t_high])
    c = np.array([C(T,S) for T in t])

    return [t, c]



#Testing
if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=[5, 4])
    axins = ax.inset_axes( np.array([0.37, 0.37, 0.6, 0.6]))#[0.005, 0.012, 25, 140])
    
    t_low = np.linspace(0.001,0.006,50)
    t_peak = np.linspace(0.006,0.02,50)
    t_high = np.linspace(0.02,0.1,50)
    t_low, t_peak, t_high = _set_temperatures(axins = axins)

    c_low = np.array([C(T,system.S) for T in t_low])

    c_peak = np.array([C(T,system.S) for T in t_peak])

    plt.title('Heat Capacity Example')
    plt.ylabel(r'$C_V$')
    plt.xlabel(r'$T$')



    c_high = np.array([C(T,system.S) for T in t_high])


    c_peak = np.array([C(T,system.S) for T in t_peak])

    

    ax.plot(np.concatenate((t_low,t_peak,t_high)), np.concatenate((c_low,c_peak,c_high)))

    # inset axes....
    axins.plot(np.concatenate((t_low,t_peak,t_high)),np.concatenate((c_low,c_peak,c_high)))
    # sub region of the original image
    x1, x2, y1, y2 = 0.002, 0.009, 5, 27
    #axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    # axins.set_xticklabels('')
    # axins.set_yticklabels('')

    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.show()







