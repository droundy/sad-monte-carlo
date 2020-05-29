import argparse, sys, glob, cbor
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


parser = argparse.ArgumentParser()
parser.add_argument('moviedir')
args = parser.parse_args()
moviedir = args.moviedir

frames = sorted(glob.glob(moviedir+'/*.cbor'))

#how to set the row and column from the cbor file?
#row = sys.argv[3]
#col = sys.argv[4]

"""
row = 0
for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        row = int(i + 1)
        break
print(row)
col = int(len(energy_data)/row)
print(col)
"""

for frame in frames:
    print('frame is', frame)
    with open(frame,'rb') as stream:
        data_loaded = cbor.load(stream)
    hist = np.array(data_loaded['bins']['histogram'])
    print('system', data_loaded['system'].keys())
    print('system E', data_loaded['system']['E'])
    print('bins', data_loaded['bins'].keys())
    print(data_loaded.keys())

    N_sites = data_loaded['system']['L']**3

    Sexcess = np.array(data_loaded['bins']['lnw'])

    NE = data_loaded['bins']['num_E']
    Nmax = data_loaded['bins']['max_N']
    moves = data_loaded['moves']

    E = np.flip(-np.arange(0, NE, 1.0))
    N = np.arange(0, Nmax+1, 1.0)
    N, E = np.meshgrid(N, E)

    Eedges = np.flip(-np.arange(-0.5, NE, 1.0))
    Nedges = np.arange(-0.5, Nmax+1, 1.0)
    Nedges, Eedges = np.meshgrid(Nedges, Eedges)

    hist = np.reshape(hist, (NE, Nmax+1))

    if hist.max() == 0:
        print('there is no data in histogram!')
        continue
    print('total data in hist is:', hist.sum())


    Sexcess = np.reshape(Sexcess, (NE, Nmax+1))
    Sexcess = Sexcess - Sexcess[-1,0]
    Sexcess[hist==0] = np.nan


    print(data_loaded['system']['E'])

    plt.figure('hist')
    plt.clf()
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    hist_to_plot = hist*1.0
    hist_to_plot[hist==0] = np.nan
    plt.pcolormesh(Nedges, Eedges, hist_to_plot)
    plt.colorbar().set_label('histogram')
    plt.title('%.3g moves' % moves)
    plt.tight_layout()
    plt.pause(1e-9)

    plt.figure('Sexcess')
    plt.clf()
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.pcolormesh(Nedges, Eedges, Sexcess)
    plt.colorbar().set_label('S_excess')
    plt.title('%.3g moves' % moves)
    plt.tight_layout()
    plt.pause(1e-9)

    Sideal = N*(1 + np.log(N_sites/N))
    S = Sexcess + Sideal

    plt.figure('entropy')
    plt.clf()
    plt.pcolor(N,E,S)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.pcolormesh(Nedges, Eedges, S)
    plt.colorbar().set_label('entropy')
    plt.tight_layout()
    plt.pause(1e-9)

    # T_inv[:-1,:-1] = (S[1:,:-1] - S[:-1,:-1])
    # T[:-1,:-1] = 1/T_inv[:-1,:-1]
    #
    # plt.figure('temperature')
    # plt.clf()
    # plt.pcolor(N,E,T, norm=LogNorm(vmin=0.01, vmax=10))
    # plt.xlabel('$N$')
    # plt.ylabel('$E$')
    # plt.pcolormesh(Nedges, Eedges, hist_to_plot)
    # plt.colorbar().set_label('temperature')
    # plt.tight_layout()
    #
    # averaged_T[:,:-1] = 2/(T_inv[:,1:] + T_inv[:,:-1])
    #
    # #what is chemical potential?
    # #chem_potential[:,:-1] = -averaged_T[:,:-1] * (S[:,1:] - S[:,:-1]) / (number_data[:,1:] - number_data[:,:-1])
    # chem_potential = np.array(data_loaded['bins']['histogram'])
    # """
    # plt.figure('excess chemical potential')
    # plt.clf()
    # #plt.title(f'{moves[t]} moves')
    # plt.pcolor(N,E,chem_potential)
    # plt.xlabel('$N$')
    # plt.ylabel('$E$')
    # plt.pcolormesh(Nedges, Eedges, hist_to_plot)
    # plt.colorbar().set_label('excess chemical potential')
    # plt.tight_layout()
    # """
    # p_exc = (T * Sexcess + chem_potential * number_data - energy_data)/N_sites**2
    # p_ideal = T * number_data / N_sites
    # pressure = p_ideal + p_exc
    #
    # plt.figure('pressure')
    # plt.clf()
    # plt.pcolor(N,E,pressure)
    # plt.xlabel('$N$')
    # plt.ylabel('$E$')
    # plt.pcolormesh(Nedges, Eedges, hist_to_plot)
    # plt.colorbar().set_label('pressure')
    # plt.tight_layout()
    #
    # #what is gibbs?
    # #gibbs_free[:,:-1] = chem_potential[:,:-1]*number_data[:,1:]
    # gibbs_free = np.zeros((col, row))
    # plt.figure('gibbs')
    # plt.clf()
    # plt.pcolor(N,E,gibbs_free)
    # plt.xlabel('$N$')
    # plt.ylabel('$E$')
    # plt.pcolormesh(Nedges, Eedges, hist_to_plot)
    # plt.colorbar().set_label('gibbs')
    # plt.tight_layout()

    plt.pause(1)





print("...and that's all, folks!")

plt.ioff()
plt.show()
