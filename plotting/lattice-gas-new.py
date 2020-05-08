import argparse, sys, glob, cbor
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()
file = args.file

frames = sorted(glob.glob(file+'/*.cbor'))


for frame in frames:
    print('frame is', frame)
    with open(frame,'rb') as stream:
        data_loaded = cbor.load(stream)
    hist = np.array(data_loaded['bins']['histogram'])
    print('system', data_loaded['system'].keys())
    print('system E', data_loaded['system']['E'])
    print('bins', data_loaded['bins'].keys())
    print(data_loaded.keys())

    NE = data_loaded['bins']['num_E']
    Nmax = data_loaded['bins']['max_N']

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
    print(Nedges.shape, Eedges.shape, hist.shape, N.shape, E.shape)
    plt.pcolormesh(Nedges, Eedges, hist)

    plt.pause(1)





print("...and that's all, folks!")

plt.ioff()
plt.show()
