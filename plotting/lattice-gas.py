import argparse, sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()
file = args.file

with open(file,'rb') as stream:
    try:
        if 'cbor' in file:
            import cbor
            data_loaded = cbor.load(stream)
        else:
            import yaml
            data_loaded = yaml.full_load(stream)
    except IOError:
        print('An error occurred trying to read the file.')


time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
hist_data = data_loaded['movies']['histogram']
hist = data_loaded['bins']['histogram']
moves = data_loaded['movies']['time']
N_sites = len(data_loaded['system']['N'])

print(moves)

number_data = np.array(data_loaded['movies']['number'])
energy_data = np.array(data_loaded['movies']['energy'])
print(energy_data)

row = 0
for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        row = int(i + 1)
        break
print(row)
col = int(len(energy_data)/row)

nlist = len(energy_data)

k = 1

energy_data.resize(col, row)
number_data.resize(col, row)

E = np.zeros((col + 1, row + 1))
E[:-1,:-1] = energy_data
N = np.zeros((col + 1, row + 1))
N[:-1,:-1] = number_data

dE = abs(E[0,0] - E[1,0])
E -= dE/2

E[-1,:] = E[-2,:] + dE
N[-1,:] = N[-2,:]

E[:,-1] = E[:,-2]
N[:,-1] = N[:,-2] + 1

N -= 0.5

S_ideal = number_data*k*(1 + np.log(N_sites/number_data))

S_ideal[number_data==0] = 0

T = np.zeros((col, row))
T_inv = np.zeros((col, row))

pressure = np.zeros((col, row))
p_ideal = np.zeros((col, row))
p_exc = np.zeros((col, row))
gibbs_free= np.zeros((col, row))

for t in range(len(entropy_data)):

    S_excess = np.array(entropy_data[t])
    S_excess.resize(col, row)
    S0_excess = S_excess[-1,0]
    S_excess = S_excess - S0_excess
    S = S_excess + S_ideal
    hist = np.array(hist_data[t])*1.0 # the multiplication converts it to floating point values
    if hist.max() == 0:
        continue
    hist.resize(col, row)
    S[hist==0] = np.nan
    S_excess[hist==0] = np.nan
    hist[hist==0] = np.nan

    plt.figure('entropy')
    plt.clf()
#   plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    plt.figure('excess entropy')
    plt.clf()
#    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S_excess)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()
    plt.figure('temperature')
    plt.clf()
#    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    plt.figure('histogram')
    plt.clf()
#    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,hist)
    plt.colorbar()

    # for i in np.arange(0,col-1,1):
    #     for j in np.arange(0,row,1):
    #         T[i][j] = dE / (S[i+1][j] - S[i][j])

    # for i in np.arange(0,col-1,1):
    #     for j in np.arange(0,row-1,1):
    #         T_inv[i][j] = (S[i+1][j] - S[i][j])/dE
    T_inv[:-1,:-1] = (S[1:,:-1] - S[:-1,:-1])/dE
    T[:-1,:-1] = 1/T_inv[:-1,:-1]

    plt.figure('temperature')
    plt.clf()
    plt.pcolor(N,E,T, norm=LogNorm(vmin=0.01, vmax=10)) #, cmap='PuBu_r'vmax=0.5, vmin=0)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    '''
    plt.figure('inverse temperature')
    plt.clf()
    plt.pcolor(N,E,T_inv, norm=LogNorm(vmin=0.01, vmax=10)) #, cmap='PuBu_r'vmax=0.5, vmin=0)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()
    '''

    """
    plt.figure('log temperature')
    plt.clf()
    plt.pcolor_log(N,E,T, norm = LogNorm(vmin = 0, vmax = 0.5))
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()
    """

    averaged_T = np.zeros_like(T)
    chem_pot = np.zeros((col, row))
    # for i in np.arange(0,col-1,1):
    #     for j in np.arange(0,row-1,1):
    #         averaged_T[i][j] = 2/(T_inv[i][j+1] + T_inv[i][j])
    #         chem_pot[i][j] = -averaged_T[i][j] * ((S[i][j+1] - S[i][j]) / (N[i][j+1] - N[i][j]))
    averaged_T[:,:-1] = 2/(T_inv[:,1:] + T_inv[:,:-1])
    chem_pot[:,:-1] = -averaged_T[:,:-1] * (S[:,1:] - S[:,:-1]) / (number_data[:,1:] - number_data[:,:-1])

    chem_pot[chem_pot==0] = np.nan
    chem_pot[T<0] = np.nan
    chem_pot[T>1] = np.nan
    chem_pot[averaged_T>1] = np.nan
    chem_pot[averaged_T<0] = np.nan

    plt.figure('excess chemical potential')
    plt.clf()
    #plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,chem_pot)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    """
    for i in np.arange(0,col-1,1):
        for j in np.arange(0,row-1,1):
            # G = chem_pot*number_data
            gibbs_free[i][j] = energy_data[i][j] - 2* T[i][j] * s_excess[i][j]

    """
    # mu = mu_ideal + mu_exc = ~?~ kT ln(N/A) + mu_exc
    #mu_ideal = T[i][j] * np.log(number_data[][]/N_sites)
    # U_exc = T*S_exc - p_exc*A + mu_exc*N # in two dimensions volume -> area


    # for i in np.arange(0,col-1,1):
    #     for j in np.arange(0,row-1,1):
    #         p_exc[i][j] = (T[i][j] * S_excess[i][j] + chem_pot[i][j] * N[i][j] - energy_data[i][j])/N_sites**2
    #         p_ideal[i,j] = T[i][j] * number_data[i][j] / N_sites
    print(T.shape, S_excess.shape, chem_pot.shape, N.shape, energy_data.shape)
    print('N is', N[0,:10], '...')
    print('number_data is', number_data[0,:10], '...')
    p_exc = (T * S_excess + chem_pot * number_data - energy_data)/N_sites**2
    p_ideal = T * number_data / N_sites

    # p = p_ideal + p_exc = kT*N/A (A = number of lattices, N = number of particles) + p_exc
    pressure = p_ideal + p_exc

    # print(pressure)

    """
    pressure[pressure==0] = np.nan
    pressure[T<0] = np.nan
    pressure[T>1] = np.nan
    pressure[averaged_T>1] = np.nan
    pressure[averaged_T<0] = np.nan
    """

    plt.figure('pressure')
    plt.clf()
    #pressure = T(ds/dv)U,N
    plt.pcolor(N,E,pressure)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    """
    plt.figure('gibbs')
    plt.clf()
    plt.pcolor(N,E,gibbs_free)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()
    """
    print('frame', t, '/', len(entropy_data))
    plt.pause(1)

print("...and that's all, folks!")

plt.ioff()
plt.show()
