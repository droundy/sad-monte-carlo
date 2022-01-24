import numpy as np
import matplotlib.pyplot as plt
import heat_capacity
import system
import glob
import os

fig, ax = plt.subplots(figsize=[5, 4], num='latest-heat-capacity')
axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])

heat_capacity.plot(system.S, ax=ax, axins=axins)
ax.indicate_inset_zoom(axins, edgecolor="black")

T_low, T_peak, T_high = heat_capacity._set_temperatures()



data_paths = sorted(glob.glob(os.path.join('tem-test', '*.npz')))
for d in data_paths:
    data = np.load(d)

    #print([(data['T'][i], c) for i, c in enumerate(data['Cv'])])

    plt.plot(data['T'], data['Cv'],'x-')

    plt.title(r'Heat Capacity for Parallel Tempering')
    plt.ylabel(r'$C_v$')
    plt.xlabel(r'$T$')

    T_mask = [(min(T_peak) <= t and max(T_high) >= t) for t in data['T']]

    C_peak = data['Cv'][T_mask]
    T_peak = data['T'][T_mask]

    axins.plot(data['T'], data['Cv'],'x-')

    plt.draw_if_interactive()
    plt.pause(1.01)