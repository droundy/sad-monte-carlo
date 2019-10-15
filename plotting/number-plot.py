#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt

def latex_float(x):
    exp = np.log10(x*1.0)
    if abs(exp) > 2:
        x /= 10.0**exp
        if ('%g' % x) == '1':
            return r'10^{%.0f}' % (exp)
        return r'%g\times 10^{%.0f}' % (x, exp)
    else:
        return '%g' % x

allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']))

my_histogram = {}
current_histogram = {}
my_entropy = {}
my_volume = {}
current_free_energy = {}
current_total_energy = {}
my_temperature = {}
my_time = {}
my_color = {}
max_iter = 0
my_gamma = {}
my_gamma_t = {}
Smin = None
fnames = sys.argv[1:]
for fname in fnames:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    my_temperature[fname] = data['T']
    current_histogram[fname] = np.array(data['bins']['histogram'])
    current_free_energy[fname] = -my_temperature[fname]*np.array(data['bins']['lnw'])
    current_free_energy[fname] = current_free_energy[fname]-current_free_energy[fname][0]
    my_volume[fname] = float(data['system']['cell']['box_diagonal']['x'])**3
    current_total_energy[fname] = np.array(data['bins']['total_energy'])
    my_color[fname] = allcolors.pop()
    my_time[fname] = np.array(data['movies']['time'])
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    my_entropy[fname] = np.array(data['movies']['lnw'])
    my_histogram[fname] = np.array(data['movies']['histogram'])
    my_gamma[fname] = np.array(data['movies']['gamma'], dtype=float)
    my_gamma_t[fname] = np.array(data['movies']['gamma_time'])
    if 'Sad' in data['method']:
        minT = data['method']['Sad']['min_T']
    if Smin is None:
        Sbest = my_entropy[fname][-1,:]
        Smin = Sbest[Sbest!=0].min() - Sbest.max()

plt.figure('gamma')
for fname in fnames:
        plt.loglog(my_gamma_t[fname], my_gamma[fname], color=my_color[fname], label=fname)
plt.legend(loc='best')
plt.xlabel('$t$')
plt.ylabel(r'$\gamma$')
# plt.ylim(1e-12, 1.1)

plt.figure('histograms')
for fname in fnames:
        plt.plot(current_histogram[fname],
                   color=my_color[fname], label=fname)
        #print(my_histogram[fname])
plt.legend(loc='best')

plt.figure('excess free energy')
for fname in fnames:
        plt.plot(current_free_energy[fname],
                   color=my_color[fname], label=fname)

plt.legend(loc='best')

plt.figure('excess internal energy')
for fname in fnames:
        plt.plot(current_total_energy[fname]/current_histogram[fname],
                   color=my_color[fname], label=fname)

plt.legend(loc='best')

plt.figure('excess entropy')
for fname in fnames:
        U = current_total_energy[fname]/current_histogram[fname]
        F = current_free_energy[fname]
        T = my_temperature[fname]
        S = (U-F)/T
        S = S-S[0]
        plt.plot(S,
                   color=my_color[fname], label=fname)

plt.figure('excess entropy/N')
for fname in fnames:
        U = current_total_energy[fname]/current_histogram[fname]
        F = current_free_energy[fname]
        T = my_temperature[fname]
        S = (U-F)/T
        S = S-S[0]
        SN = np.arange(0, len(S), 1)
        plt.plot((np.pi/6)*SN/my_volume[fname],S/SN,
                   color=my_color[fname], label=fname)
        plt.xlabel(r'$\eta$')


plt.legend(loc='best')

plt.figure('excess internal energy/N')
for fname in fnames:
        U = current_total_energy[fname]/current_histogram[fname]
        UN = np.arange(0, len(U), 1)
        plt.plot((np.pi/6)*UN/my_volume[fname],U/UN,
                   color=my_color[fname], label=fname)

plt.figure('Pressure')
for fname in fnames:
        U = current_total_energy[fname]/current_histogram[fname]
        F = current_free_energy[fname]
        V = my_volume[fname]
        T = my_temperature[fname]
        N = len(F)
        p = np.zeros(N-1)
        p_exc = np.zeros(N-1)
        for i in range(0,N-1):
                u = F[i+1]-F[i] # dN = 1
                p_exc[i] = (-F[i]+u*(i+.5))/V
                p[i] = (-F[i]+u*(i+.5))/V+(i+.5)*T/V
        UN = np.arange(0.5, N-1, 1)
        #print(len(UN), len(p))
        plt.ylabel('Pressure')
        plt.xlabel(r'$\eta$')
        plt.plot((np.pi/6)*UN/my_volume[fname],p,
                   color=my_color[fname], label=fname)
        plt.plot((np.pi/6)*UN/my_volume[fname],p_exc,'--',
                   color=my_color[fname], label=fname + ' pexc')


plt.figure('Gibbs')
for fname in fnames:
        U = current_total_energy[fname]/current_histogram[fname]
        F = current_free_energy[fname]
        V = my_volume[fname]
        T = my_temperature[fname]
        N = len(F)
        p = np.zeros(N-1)
        p_exc = np.zeros(N-1)
        for i in range(0,N-1):
                u = F[i+1]-F[i] # dN = 1
                p_exc[i] = (-F[i]+u*(i+.5))/V
                p[i] = (-F[i]+u*(i+.5))/V+(i+.5)*T/V
        G = np.zeros(N-2)
        p_integer = np.zeros(N-2)
        for j in range(1,N-2):
                p_integer[j] = (p[j]+p[j+1])/2
                G[j] = F[j] + V*p_integer[j]
        plt.ylabel('Gibbs')
        plt.xlabel('Pressure')
        plt.plot(p_integer,G,
                   color=my_color[fname], label=fname)
plt.legend(loc='best')

all_mu = np.linspace(-10, 2, 30000)
# nQ = (mkT/2pi hbar^2)^1.5
nQ = 0.001 # HOKEY

plt.figure('Grand Uexc')
for fname in fnames:
        Uexc_N = current_total_energy[fname]/current_histogram[fname]
        Fexc_N = current_free_energy[fname]
        V = my_volume[fname]
        T = my_temperature[fname]
        beta = 1/T
        Nmax = len(Fexc_N)-1
        N_N = np.arange(0, Nmax+1, 1)
        # print(len(Fexc_NR))
        Fid_N = N_N*T*np.log(N_N/V/nQ) - N_N*T
        Fid_N[0] = 0
        # print(len(Fid_N))
# # you need to figure out the indexing for Z grand. Make sure it sums over all number. Get a 60 vs 61 index error.
        Grand_Uexc = np.zeros_like(all_mu)
        Grand_N = np.zeros_like(all_mu)
        for i in range(len(all_mu)):
            mu = all_mu[i]
            # Zgrand = \sum_N e^{-\beta(Fid(N) + Fexc_N - mu N)}
            Zgrand_exponents = -beta*(Fid_N+Fexc_N-mu*N_N)
            # print(Zgrand_exponents)
            offset = Zgrand_exponents.max()
            Zgrand_exponents -= offset
            Zgrand = np.exp(Zgrand_exponents).sum()
            Grand_Uexc[i] = (Uexc_N*np.exp(Zgrand_exponents)).sum()/Zgrand
            Grand_N[i] = (N_N*np.exp(Zgrand_exponents)).sum()/Zgrand
        #U_ideal = C_V*T*N*k : where C_V i think is 3
        C_V = 3
        Grand_Uideal = C_V*T*Grand_N
        plt.plot(Grand_N, Grand_Uexc + Grand_Uideal,'-',
                   color=my_color[fname], label=fname)
        plt.plot(current_total_energy[fname]/current_histogram[fname], ':',
                   color=my_color[fname], label='canonical '+fname)
        plt.ylabel('Uexc')
        plt.xlabel('N')
        plt.legend(loc='best')

plt.figure('Grand S_Excess')
for fname in fnames:
        Uexc_N = current_total_energy[fname]/current_histogram[fname]
        Fexc_N = current_free_energy[fname]
        T = my_temperature[fname]
        Sexc_N = (U-F)/T
        Sexc_N = Sexc_N-Sexc_N[0]
        V = my_volume[fname]
        beta = 1/T
        Nmax = len(Fexc_N)-1
        N_N = np.arange(0, Nmax+1, 1)
        Grand_Sexc = np.zeros_like(all_mu)
        Grand_N = np.zeros_like(all_mu)
        Grand_Sideal = np.zeros_like(all_mu)
        for i in range(len(all_mu)):
            mu = all_mu[i]
            # Zgrand = \sum_N e^{-\beta(Fid(N) + Fexc_N - mu N)}
            Zgrand_exponents = -beta*(Fid_N+Fexc_N-mu*N_N)
            offset = Zgrand_exponents.max()
            Zgrand_exponents -= offset
            Zgrand = np.exp(Zgrand_exponents).sum()
            Grand_Sexc[i] = (Sexc_N*np.exp(Zgrand_exponents)).sum()/Zgrand
            Grand_Sideal[i] = (Grand_Uideal[i]+mu*Grand_N[i])
            Grand_N[i] = (N_N*np.exp(Zgrand_exponents)).sum()/Zgrand
        # Grand_Sideal = ????
        plt.plot(Grand_N, Grand_Sexc,':',
                    color=my_color[fname], label=fname)
        plt.plot(Grand_N, Grand_Sideal,':',
                    color=my_color[fname], label='ideal')
        plt.ylabel('Sexc')
        plt.xlabel('N')

# print(len(Grand_Uexc))

plt.figure('Grand P_exc')
for fname in fnames:
    mu = all_mu[i]
    Uexc_N = current_total_energy[fname]/current_histogram[fname]
    Fexc_N = current_free_energy[fname]
    T = my_temperature[fname]
    V = my_volume[fname]
    Grand_Pexc = np.zeros_like(all_mu)
    for i in range(len(all_mu)):
        Grand_Pexc[i] = (T*Grand_Sexc[i] + mu * Grand_N[i] - Grand_Uexc[i])/V
    plt.plot(Grand_N, Grand_Pexc,':',
                color=my_color[fname], label=fname)
    p_ideal = T*Grand_N/V
    plt.plot(Grand_N, p_ideal,'-',
                color=my_color[fname], label='ideal')

plt.figure('Grand G vs P')
for fnamme in fnames:
    mu = all_mu[i]
    Grand_G = np.zeros_like(all_mu)
    # Grand_G = mu*Grand_N
    for i in range(len(all_mu)):
        Grand_G[i] = mu*Grand_N[i]
    p_ideal = T*Grand_N/V
    plt.plot(Grand_Pexc + p_ideal, Grand_G,'.:',
                color=my_color[fname], label=fname)

plt.xlabel('$p$')
plt.ylabel('$G$')

plt.show()
