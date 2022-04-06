import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import heat_capacity
from styles import linestyle
import system
import glob
import os
import find_phase_transition

def normalize_S(S, E=None):
    if E is None:
        E = np.linspace(-system.h_small, 0, 1000)

    dE = np.diff(E)
    E = 0.5*(E[1:] + E[:-1])
    S = S - max(S)
    total = np.sum(0.5*(np.exp(S)[1:] + np.exp(S)[:-1])*dE)
    return S - np.log(total)

def make_plots(Ts, mean_E, mean_Esqr, label):
    if label == 'exact':
        linestyle = ':'
    else:
        linestyle = '--'
    refined_Ts = np.linspace(Ts[np.argmin(abs(Ts-0.005))], Ts[np.argmin(abs(Ts-0.25))], 30*len(Ts))
    mean_E_func = interp1d(Ts, mean_E, fill_value='extrapolate', kind='cubic')
    mean_E_sqr_func = interp1d(Ts, mean_Esqr, fill_value='extrapolate', kind='cubic')
    variance = mean_E_sqr_func(refined_Ts) - mean_E_func(refined_Ts)**2
    #var_func = interp1d(Ts, variance, fill_value='extrapolate', kind='cubic')
    
    C = (mean_E_sqr_func(refined_Ts) - mean_E_func(refined_Ts)**2)/ (refined_Ts**2)
    plt.figure('variance')
    plt.vlines(find_phase_transition.actual_T,min(variance),max(variance), alpha=0.25)

    plt.plot(refined_Ts, variance, label=label, linestyle=linestyle)
    plt.ylabel(r'Var$(E)$')
    plt.xlabel(r'$T$')
    plt.title('Variance of Energy')
    plt.xlim([0, 0.08])


    plt.legend()

    plt.figure('mean_E')
    plt.vlines(find_phase_transition.actual_T,min(mean_E),max(mean_E), alpha=0.25)
    
    plt.plot(Ts, mean_E, label=label, linestyle=linestyle)
    plt.ylabel(r'$\langle E \rangle _T$')
    plt.xlabel(r'$T$')
    plt.title('Mean of Energy')
    plt.xlim([0, 0.08])


    plt.legend()

    plt.figure('mean_Esqr')
    plt.vlines(find_phase_transition.actual_T,min(mean_Esqr),max(mean_Esqr), alpha=0.25)
    

    plt.plot(Ts, mean_Esqr, label=label, linestyle=linestyle)
    plt.ylabel(r'$\langle E^2 \rangle _T$')
    plt.xlabel(r'$T$')
    plt.title('Mean of the Squared Energy')
    plt.xlim([0, 0.08])


    plt.legend()

    plt.figure('C_V')
    plt.vlines(find_phase_transition.actual_T,min(C),max(C), alpha=0.25)

    plt.plot(refined_Ts, C, label=label, linestyle=linestyle)
    plt.ylabel(r'$C_V(T)$')
    plt.xlabel(r'$T$')
    plt.title('Heat Capacity')
    plt.xlim([0, 0.08])
    plt.legend()
    


T_low,T_peak,T_high = heat_capacity._set_temperatures()
movie = False
if movie:
    fig, ax = plt.subplots(figsize=[5, 4], num='latest-heat-capacity')
    axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])

    heat_capacity.plot(system.S, ax=ax, axins=axins)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    data_paths = sorted(glob.glob(os.path.join('tem+T-trans-1+barrier-0+seed-1', '*.npz')))
    for d in data_paths:
        
        data = np.load(d)

        data_10 = np.load(d)

        variance_10 = interp1d(data_10['T'], data_10['var_E'], fill_value='extrapolate', kind='cubic')
        cv_10 = lambda t: variance_10(t)/(t**2)

        Ts = np.linspace(min(data_10['T']),max(data_10['T']),1000)

        plt.plot(Ts, cv_10(Ts), label = 10)
        axins.plot(Ts, cv_10(Ts))

        plt.xlim(0,0.25)

        #plt.plot(data['T'], data['Cv'],'x-')
        #plt.plot(data['T'], data['Cv_old'], 'o-')

        plt.title(r'Heat Capacity for Parallel Tempering')
        plt.ylabel(r'$C_v$')
        plt.xlabel(r'$T$')

        T_mask = [(min(T_peak) <= t and max(T_high) >= t) for t in data['T']]

        C_peak = data['Cv'][T_mask]
        T_peak = data['T'][T_mask]

        #axins.plot(data['T'], data['Cv'],'x-')

        plt.draw_if_interactive()
        plt.pause(1.01)
    plt.show()
else:
    data_paths_50 = sorted(glob.glob(os.path.join('tem+T-trans-1+barrier-0+seed-1', '*.npz')))
    data_50 = np.load(data_paths_50[-1])
    for k in data_50.keys():
        print(k)


    variance_50 = interp1d(data_50['T'], data_50['var_E'], fill_value='extrapolate', kind='linear')

    cv_50 = lambda t: variance_50(t)/(t**2)

    T_low,T_peak,T_high = heat_capacity._set_temperatures(Tmax = .5)
    Ts = np.concatenate((T_low[:-1],T_peak,T_high[1:]))

    mean_E_exact=[]
    mean_Esqr_exact=[]
    for t in Ts:
        C, mean_E, mean_Esqr = heat_capacity.C_E_Esqrd(t, system.S)
        mean_E_exact.append(mean_E)
        mean_Esqr_exact.append(mean_Esqr)

    mean_E_exact=np.array(mean_E_exact)
    mean_Esqr_exact=np.array(mean_Esqr_exact)
    variance_exact = mean_Esqr_exact - mean_E_exact**2

    var_exact = interp1d(Ts, variance_exact, fill_value='extrapolate', kind='linear')
    mean_E_exact_f = interp1d(Ts, mean_E_exact, fill_value='extrapolate', kind='linear')
    C_exact = lambda t: var_exact(t)/t**2

    plt.figure('mean-E')
    plt.plot(Ts, mean_E_exact, 'x-')

    var_tempering = interp1d(data_50['T'], data_50['var_E'], fill_value='extrapolate', kind='linear')
    mean_E_tempering_f = interp1d(data_50['T'], data_50['mean_E'], fill_value='extrapolate', kind='linear')
    C_tempering = lambda t: var_tempering(t)/t**2

    T_int = np.linspace(0.0001, .5, 10000)[1:]
    dT = T_int[1] - T_int[0]

    #C_int_exact = C_exact(T_int)*dT/T_int
    #C_int_tempering = C_tempering(T_int)*dT/T_int
    plt.figure('C')
    plt.plot(T_int, C_exact(T_int))

    #S_int = np.flip(np.cumsum(np.flip(-C_int_exact)))
    #S_int_2 = np.cumsum(C_int_exact)
    #S_tempering = np.flip(np.cumsum(np.flip(-C_int_tempering)))
    #plt.figure('S')

    #E = np.linspace(-system.h_small, 0, 10000)
    #plt.plot(mean_E_exact_f(T_int), normalize_S(S_int_2, mean_E_exact_f(T_int)) , label= 'exact by integrating Cv')
    #plt.plot(E, normalize_S(system.S(E), E), ':', label='exact')
    #plt.plot(mean_E_tempering_f(T_int), S_tempering, label = 'Par Temp')
    #lt.title(r'$S(E)$ from integrating heat capacity')
    #plt.legend()

    make_plots(Ts,mean_E_exact,mean_Esqr_exact,'exact')
    make_plots(data_50['T'],data_50['mean_E'], data_50['var_E'] + data_50['mean_E']**2, 'Paralell Tempering')



plt.show()