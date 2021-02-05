import numpy as np
import scipy.optimize as optimize
import yaml

def read_file(base):
    energy_boundaries = np.loadtxt(base+'-energy-boundaries.dat')
    mean_energy = np.loadtxt(base+'-mean-energy.dat')
    lnw = np.loadtxt(base+'-lnw.dat')
    with open(base+'-system.dat') as f:
        system = yaml.safe_load(f)

    if energy_boundaries[0] < energy_boundaries[-1]:
        energy_boundaries = np.flip(energy_boundaries)
        mean_energy = np.flip(mean_energy)
        lnw = np.flip(lnw)
    lnw -= lnw.max()
    lnw -= np.log(np.sum(np.exp(lnw)))
    return energy_boundaries, mean_energy, lnw, system

def step_entropy(energy_boundaries, mean_energy, lnw):
    step_entropy = []
    step_energy = []
    for i in range(len(energy_boundaries)-1):
        step_energy.append(energy_boundaries[i])
        step_energy.append(energy_boundaries[i+1])
        Shere = lnw[i+1] - np.log(energy_boundaries[i] - energy_boundaries[i+1])
        step_entropy.append(Shere)
        step_entropy.append(Shere)
    step_energy = np.flip(np.array(step_energy))
    step_entropy = np.flip(np.array(step_entropy))
    def entropy(E):
        if len(step_energy) < 3:
            return np.zeros_like(E)
        # return np.interp(E, step_energy, step_entropy, left=step_entropy[0], right=step_entropy[-1])
        return np.interp(E, step_energy, step_entropy, left=-30, right=-30)
    return entropy, 1*step_energy, 1*step_entropy

def fn_for_beta(x, meanE_over_deltaE):
    if x < 1e-14:
        return 0.5*x - x*meanE_over_deltaE
    return x/(1-np.exp(-x)) - 1 - x*meanE_over_deltaE
    
def find_beta_deltaE(meanE_over_deltaE):
    # x = np.linspace(-100,100,10000)
    # plt.plot(x, np.vectorize(fn_for_beta)(x, meanE_over_deltaE))
    # plt.show()
    if meanE_over_deltaE == 0:
        x0 = 1e-6
        x1 = -1e-6
    else:
        x0 = 2*meanE_over_deltaE
        x1 = meanE_over_deltaE
    sol = optimize.root_scalar(fn_for_beta, args=(meanE_over_deltaE), x0 = x0, x1 = x1)
    # print(sol)
    return sol.root

def find_entropy_from_beta_and_lnw(beta, lnw, deltaE):
    if abs(beta*deltaE) < 1e-14:
        return lnw - np.log(deltaE)
    return lnw - np.log(deltaE) - np.log((np.exp(beta*deltaE)-1)/(beta*deltaE))

def linear_entropy(energy_boundaries, mean_energy, lnw):
    step_entropy = []
    step_energy = []
    total_energy_delta = 2*(energy_boundaries[0] - energy_boundaries[-1])
    # First we include the unbounded high-energy case
    # step_energy.append(energy_boundaries[0] + total_energy_delta)
    # step_energy.append(energy_boundaries[0])
    sigup = (mean_energy[0] - energy_boundaries[0])*np.sqrt(np.pi/2)
    Sup = lnw[0] - np.log(sigup*np.sqrt(np.pi/2))
    # step_entropy.append(S0-total_energy_delta/kT)
    # step_entropy.append(S0)

    # Now consider all the middle cases
    for i in range(len(energy_boundaries)-1):
        step_energy.append(energy_boundaries[i])
        step_energy.append(energy_boundaries[i+1])

        deltaE = energy_boundaries[i] - energy_boundaries[i+1]
        meanE = mean_energy[i+1]
        beta = find_beta_deltaE((meanE - energy_boundaries[i+1])/deltaE)/deltaE
        S0 = find_entropy_from_beta_and_lnw(beta, lnw[i+1], deltaE)
        if np.isnan(beta+S0):
            step_entropy.append(lnw[i+1]-np.log(deltaE))
            step_entropy.append(lnw[i+1]-np.log(deltaE))
        else:
            step_entropy.append(S0+beta*deltaE)
            step_entropy.append(S0)
    # Now let's do the low-energy unbounded bin
    kT = energy_boundaries[-1] - mean_energy[-1]
    S0 = lnw[-1] - np.log(kT)
    if not np.isnan(kT) and not np.isnan(S0):
        step_entropy.append(S0)
        step_entropy.append(S0-total_energy_delta/kT)
        step_energy.append(energy_boundaries[-1])
        step_energy.append(energy_boundaries[-1] - total_energy_delta)

    step_energy = np.flip(step_energy)
    step_entropy = np.flip(step_entropy)
    def entropy(E):
        if np.isnan(Sup) or np.isnan(sigup):
            return np.interp(E, step_energy, step_entropy, left=step_entropy[0], right=-30)
        else:
            return np.interp(E, step_energy, step_entropy, left=step_entropy[0], right=0) \
                + np.heaviside(E-energy_boundaries[0], 0)*(Sup-(E-energy_boundaries[0])**2/(2*sigup**2))
    return entropy, 1*step_energy, 1*step_entropy
