import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import heat_capacity
from styles import linestyle
import system
import glob
import os

fig, ax = plt.subplots(figsize=[5, 4], num='latest-heat-capacity')
axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])

T_exact, C_exact = heat_capacity.plot(system.S, ax=ax, axins=axins)
ax.indicate_inset_zoom(axins, edgecolor="black")

T_low, T_peak, T_high = heat_capacity._set_temperatures()
movie = True
if movie:

    data_paths = sorted(glob.glob(os.path.join('tem-T-50', '*.npz')))
    for d in data_paths:
        data = np.load(d)

        data_10 = np.load(d)

        variance_10 = interp1d(data_10['T'], data_10['var_E'], fill_value='extrapolate', kind='cubic')
        cv_10 = lambda t: variance_10(t)/(t**2)

        Ts = np.linspace(min(data_10['T']),max(data_10['T']),100000)

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
    data_paths_10 = sorted(glob.glob(os.path.join('tem-T-10', '*.npz')))
    data_10 = np.load(data_paths_10[-1])
    data_paths_50 = sorted(glob.glob(os.path.join('tem-T-50', '*.npz')))
    data_50 = np.load(data_paths_50[-1])


    variance_10 = interp1d(data_10['T'], data_10['var_E'], fill_value='extrapolate', kind='cubic')
    variance_50 = interp1d(data_50['T'], data_50['var_E'], fill_value='extrapolate', kind='cubic')

    cv_10 = lambda t: variance_10(t)/(t**2)
    cv_50 = lambda t: variance_50(t)/(t**2)

    Ts = np.linspace(min(data_10['T']),max(data_10['T']),100000)

    plt.plot(Ts, cv_10(Ts), label = 10)
    axins.plot(Ts, cv_10(Ts))
    plt.plot(data_10['T'], data_10['var_E']/(data_10['T']**2), 'x')
    axins.plot(data_10['T'], data_10['var_E']/(data_10['T']**2), 'x')

    plt.plot(Ts, cv_50(Ts), label = 50)
    axins.plot(Ts, cv_50(Ts))
    plt.plot(data_50['T'], data_50['var_E']/(data_50['T']**2), 'x')
    axins.plot(data_50['T'], data_50['var_E']/(data_50['T']**2), 'x')

    plt.legend()
    plt.xlim(0,0.25)

    plt.figure('variance')

    plt.plot(T_exact, C_exact*T_exact**2, label='exact', linestyle=':')

    plt.plot(Ts, variance_10(Ts), label = 10)
    plt.plot(data_10['T'], data_10['var_E'], 'x')
    plt.plot(Ts, variance_50(Ts), label = 50)
    plt.plot(data_50['T'], data_50['var_E'], 'x')

    plt.legend()
    plt.xlim(0,0.25)

plt.show()