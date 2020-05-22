import argparse, sys, glob, cbor
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


parser = argparse.ArgumentParser()
parser.add_argument('moviedir')
args = parser.parse_args()
moviedir = args.moviedir

frames = sorted(glob.glob(moviedir+'/*.cbor'))


for frame in frames:
    print('frame is', frame)
    with open(frame,'rb') as stream:
        data_loaded = cbor.load(stream)
    hist = np.array(data_loaded['bins']['histogram'])
    print('system', data_loaded['system'].keys())
    print('system E', data_loaded['system']['E'])
    print('bins', data_loaded['bins'].keys())
    print(data_loaded.keys())

    Sexcess = np.array(data_loaded['bins']['lnw'])
    Sexcess = Sexcess - Sexcess[1,0]

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

    plt.figure('S_excess')
    plt.clf()
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.pcolormesh(Nedges, Eedges, S_excess)
    plt.colorbar().set_label('S_excess')
    plt.title('%.3g moves' % moves)
    plt.tight_layout()

    plt.pause(1)





print("...and that's all, folks!")

plt.ioff()
plt.show()
