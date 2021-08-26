import numpy as np
import matplotlib.pylab as plt
from scipy.special import erf
from scipy.special import erfi
from scipy import optimize


print('START -----------------------')

'''
Test Data
'''
# variables
E_theory = np.linspace(-5,5,1000) #E_theory is (a somewhat) continous Energy specctrum
E_exp = np.arange(-5, 6, 1) # E_exp is ten points that I'm pretending were given
sigma = 1.0 # width of the guassian
kb = 1.38064852*10**-23 # bolzman constant

# density of states function - a simple guassian function

def DOS(E):
    return 1/(np.sqrt(np.pi)*sigma)*np.exp(-E**2/(2*sigma**2))

# The solved integrals of the Density of states function give the known moments
# I'm calling the moments, I dont remember if this is the right term for them I just
# recall it started with an m. 

def known_moments(E):
    EL = E[1:]
    ER = E[:-1]
    known_weights = (erf(EL/(np.sqrt(2)*sigma))-erf(ER/np.sqrt(2)))/2
    a = (2*sigma)/np.sqrt(2*np.pi)
    b = np.exp(-ER/(2*sigma**2)) - np.exp(-EL/(2*sigma**2))
    c = erf(EL/(np.sqrt(2)*sigma)) - erf(ER/(np.sqrt(2)*sigma))
    known_E_avg = a*(b/c)
    print('known_E_avg', known_E_avg)
    d = ER*np.exp(-ER/(2*sigma**2)) - EL*np.exp(-EL/(2*sigma**2))
    known_E_sqrd_avg = (2*sigma)*(np.sqrt(2/np.pi)*(d/c**2)+(sigma/c))
    known_E_uncertainty = np.sqrt(1/known_weights*((sigma*d)/2 + (sigma**2*c/2) - (sigma**2*b/2*np.pi)))
    return [known_weights, known_E_avg, known_E_sqrd_avg, known_E_uncertainty]

print('E_exp', E_exp)
print(known_moments(E_exp))

'''
Calculate Moments
'''
# Constant Initial Guess

A_guess = 1
SL_guess = 1
SR_guess = 1
eps_guess = 1

# Notes:
# SL = Si
# SR = Si-1
# EL - Ei
# ER = Ei-1
# A must be sorted into positive or negative before applying erfi or erf

# A = x[0], SL = x[1], SR = x[2], epsilon = x[3]

def Moments(x):
    EL = 1
    # coefficients
    A = x[0]
    SR = x[1]
    SL = x[2]
    epsilon = x[3]  # EL - ER
    gamma = x[0]*epsilon**2 # A = epsilon^2
    delta_S = SL - SR
    # calculate moments
    if x[0] > 0:
        # split up function for ease of checking
        a = np.sqrt(np.pi)*np.exp(SL)
        b = np.exp(((gamma-delta_S)**2)/(gamma*4))
        c = 2*np.sqrt(gamma)
        d = erfi((delta_S + gamma)/(c))
        f = erfi((delta_S - gamma)/(c))
        g = erfi(d-f)
        weight = (a*b*g)/c
        h = 1-np.exp(-delta_S)
        # cannot express as eps but also cannot have 4 unknowns !!!
        E_avg = (epsilon/gamma)*(((delta_S + gamma)/2)-(h/g)) + EL
        k = (np.sqrt(np.pi)/(8*gamma**(5/2)))*((delta_S + gamma)**2 + 2*gamma)
        p = np.exp(((gamma+delta_S)**2)/(gamma*4))
        Theta = k*(p*d + b*f)
        q = (epsilon*np.exp(SL)/weight)
        r = 2*EL*(EL - E_avg)
        E_sqrd_avg = q*Theta + r + EL**2
        E_uncertainty = np.sqrt(E_sqrd_avg**2 - E_avg**2)
        return [weight, E_avg, E_sqrd_avg, E_uncertainty]
    else:
        print('A < 0')
        return [np.NaN, np.NaN, np.NaN, np.NaN]
        # A < 0 was removed because when A is negative np.sqrt(gamma) is imaginary !!!

'''
Root Finding
'''
def find_residuals(x, which_bin):
    # temporary values
    M = Moments(x)
    w = M[0]
    E_avg = M[1]
    E_sqrd_avg = M[2]
    E_uncertainty = M[3]
    w_matrix = w*np.ones(len(E_exp)-1)
    E_avg_matrix = w*np.ones(len(E_exp)-1)
    E_sqrd_avg_matrix = w*np.ones(len(E_exp)-1)
    E_uncertainty_matrix = w*np.ones(len(E_exp)-1)
    # know values
    K = known_moments(E_exp)
    known_w = K[0]
    known_E_avg = K[1]
    known_E_sqrd_avg = K[2]
    known_E_uncertainty = K[3]
    w_diff = w_matrix - known_w
    E_avg_diff = E_avg_matrix - known_E_avg
    E_sqrd_avg_diff = E_sqrd_avg_matrix - known_E_sqrd_avg
    E_uncertainty_diff = E_uncertainty_matrix - known_E_uncertainty
    # check values
    print('w_matrix - known_w = ', w_diff)
    print('sum = ',sum(w_diff))
    print('E_avg_matrix - known_E_avg = ', E_avg_diff)
    print('sum = ', sum(E_avg_diff))
    print('E_sqrd_avg_matrix - known_E_sqrd_avg = ', E_uncertainty_diff)
    print('sum = ', sum(E_uncertainty_diff))
    print('E_uncertainty_matrix - known_E_uncertainty = ', E_uncertainty_diff)
    print('sum = ', sum(E_uncertainty_diff))
    return w_diff[which_bin], E_avg_diff[which_bin], E_sqrd_avg_diff[which_bin], E_uncertainty_diff[which_bin]

# function for roots coppied from separte document to check functioning on simmple expression

def roots(which_bin, x0_guess, x1_guess, x2_guess,x3_guess):
    sol = optimize.root(find_residuals,[x0_guess, x1_guess, x2_guess, x3_guess], args=(which_bin))
    return sol.x

def find_coefficients(which_bin, A_guess, SL_guess, SR_guess, eps_guess):
    # use rootfinding to determine roots for weights
    A_guess, SL_guess, SR_guess, eps_guess = roots(which_bin, A_guess, SL_guess, SR_guess, eps_guess)
    # check the residuals
    w_diff, E_avg_diff, E_sqrd_avg_diff, E_uncertainty_diff = find_residuals(
        [A_guess, SL_guess, SR_guess, eps_guess], which_bin)
    print('coef = ',A_guess, SL_guess,SR_guess,eps_guess)
    return A_guess, SL_guess, SR_guess, eps_guess


print('coefficients are', find_coefficients(0, A_guess, SL_guess, SR_guess, eps_guess))

'''
Plots
'''

# # plot Density of states
# fig1 = plt.figure(1)
# fig1.set_size_inches(4, 4)
# plt.plot(E_theory, DOS(E_theory), color='black')
# plt.plot(E_exp, DOS(E_exp), 'o', color='red')
# plt.xlabel('E (J)')
# plt.ylabel('Density of States')
# plt.tight_layout()

# plot Average Energy
# fig2 = plt.figure(2)
# fig2.set_size_inches(4, 4)
# E_array = (E_theory[:-1] + E_theory[1:])/2
# SL = Entropy(E_exp)[:-1]
# SR = Entropy(E_exp)[1:]
# EL = E_exp[:-1]
# ER = E_exp[1:]
# A = 1
#print(np.shape(SL), np.shape(SR), np.shape(A), np.shape(EL), np.shape(ER))
# w1, E1 = Moments(SL, SR, A, EL, ER)
# print('E = ', E1)
# print('w = ', w1)
# plt.plot(E_exp[1:], E1, 'o', color = 'blue')
# plt.plot(E_array, Average_Energy(E_theory), color='black')
# #plt.plot(E_exp[1:], Average_Energy(E_exp), 'o', color='red')
# plt.xlabel('E (J)')
# plt.ylabel('Average Energy (J)')
# plt.tight_layout()

# # plot entropy
# fig3 = plt.figure(3)
# fig3.set_size_inches(4, 4)
# E_array = (E_theory[:-1] + E_theory[1:])/2
# plt.plot(E_theory, Entropy(E_theory), color='black')
# plt.plot(E_exp, Entropy(E_exp), 'o', color='red')
# plt.xlabel('E (J)')
# plt.ylabel('Entropy')
# plt.tight_layout()
plt.show()



print('END -------------------------')
