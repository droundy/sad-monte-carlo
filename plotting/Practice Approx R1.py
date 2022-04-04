import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special

"""setup estimation data"""

A = 3.1
B = 0.2
C = 0.001

Esize = 100
Estep = 0.1
E = np.arange(0,Esize,Estep)
dE = Estep

realS = A + B*E + C*E**2
realDos = np.exp(realS)

binnum = 10
bstep = Estep*binnum
bins = np.arange(0+bstep/2,Esize+bstep/2,bstep) #points at the centers, not the dividers


#calculating W
W = np.zeros_like(bins)
for i in range(0,len(bins)):
    sumW = 0
    for x in range(0,len(E)): # for j in range(binnum): sumW += realDos[i*binnum + j]*dE ...
        if np.abs(E[x]-bins[i]) < bstep/2:
            sumW += realDos[x] * dE
            #print(i,x)
        if np.abs(E[x]-bins[i]) == bstep/2:
            sumW += realDos[x] * dE/2
    W[i] = sumW


#calculating Ave
Ave = np.copy(bins)
E_times_Dos = E*realDos
for i in range(0,len(bins)):
    sum = 0
    for j in range(0,len(E)):
        if np.abs(E[j]-bins[i]) < bstep/2:
            sum += E[j]*realDos[j] * dE
        if np.abs(E[j]-bins[i]) == bstep/2:
            sum += E[j]*realDos[j] * dE/2
    # sum += E_times_Dos[abs(E - i) < bstep/2].sum()*dE
    Ave[i] = sum/W[i]


#calculating Ave2
Ave2 = np.copy(bins)
for i in bins:
    sum = 0
    for x in E:
        if np.abs(x-i) < bstep/2:
            sum += x*x*realDos[np.where(E==x)] * dE
        if np.abs(x-i) == bstep/2:
            sum += x*x*realDos[np.where(E==x)] * dE/2
    Ave2[np.where(bins==i)] = sum/W[np.where(bins==i)]



'''Testing Analytic Solutions:'''
# calcW = np.zeros_like(bins)
# calcAve = np.copy(bins)
# calcAve2 = np.copy(bins)


'''Note, calcW and W vary more than wanted when Dos increases steeply,
but using smaller bins showed that the inaccuracy was in the W approximation size.'''

def calcW(input_index,A,B,C):
    #Other imported info: bins[], bstep
    x = input_index
    E0 = bins[x] - bstep*0.5
    E1 = bins[x] + bstep*0.5
    if C < 0:
        F3 = B/(2*np.sqrt(-C))
        F1 = special.erf(np.sqrt(-C)*E1-F3) - special.erf(np.sqrt(-C)*E0-F3)
        calcW = 0.5*np.sqrt(-np.pi/C)*np.exp(A+F3*F3)*F1
    else:
        F3 = B/(2*np.sqrt(C))
        F1 = special.erfi(np.sqrt(C)*E1+F3) - special.erfi(np.sqrt(C)*E0+F3)
        calcW = 0.5*np.sqrt(np.pi/C)*np.exp(A-F3*F3)*F1
    return calcW

def calcAve(input_index,A,B,C):
    #Other imported info: bins[], bstep
    x = input_index
    E0 = bins[x] - bstep*0.5
    E1 = bins[x] + bstep*0.5
    if C < 0:
        F3 = B/(2*np.sqrt(-C))
        F1 = special.erf(np.sqrt(-C)*E1-F3) - special.erf(np.sqrt(-C)*E0-F3)
        comp1 = np.exp(B*E0+C*E0*E0) - np.exp(B*E1+C*E1*E1)
        calcAve = 1/(np.sqrt(-C*np.pi)*np.exp(F3**2))*comp1/F1+F3/np.sqrt(-C)
    else:
        F3 = B/(2*np.sqrt(C))
        F1 = special.erfi(np.sqrt(C)*E1+F3) - special.erfi(np.sqrt(C)*E0+F3)
        comp1 = np.exp(B*E0+C*E0*E0) - np.exp(B*E1+C*E1*E1)
        calcAve = -1/(np.sqrt(C*np.pi))*np.exp(F3**2)*comp1/F1-F3/np.sqrt(C)
    return calcAve

def calcAve2(input_index,A,B,C):
    #Other imported info: bins[], bstep
    x = input_index
    E0 = bins[x] - bstep*0.5
    E1 = bins[x] + bstep*0.5
    if C < 0:
        F3 = B/(2*np.sqrt(-C))
        F1 = special.erf(np.sqrt(-C)*E1-F3) - special.erf(np.sqrt(-C)*E0-F3)
        comp2 = np.exp(B*E0+C*E0*E0)*(B-2*C*E0) - np.exp(B*E1+C*E1*E1)*(B-2*C*E1)
        calcAve2 = 1/(2*(-C)**1.5*np.pi**0.5*np.exp(F3**2))*comp2/F1+(B**2-2*C)/(4*C**2)
    else:
        F3 = B/(2*np.sqrt(C))
        F1 = special.erfi(np.sqrt(C)*E1+F3) - special.erfi(np.sqrt(C)*E0+F3)
        comp2 = np.exp(B*E0+C*E0*E0)*(B-2*C*E0) - np.exp(B*E1+C*E1*E1)*(B-2*C*E1)
        calcAve2 = 1/(2*(C)**1.5*np.pi**0.5)*np.exp(F3**2)*comp2/F1+(B**2-2*C)/(4*C**2)
    return calcAve2



'''Rough approx estimate for A,B,C'''
#Coming soon




bindex = range(len(bins))

plt.plot(E,realDos,label="real Dos")
plt.plot(bins,W,label="W from Dos")
plt.plot(bins,calcW(bindex,A,B,C),label="calcW")
plt.legend()
plt.figure()

# plt.plot(E,realDos,label="realDos")
# plt.plot(bins,W,'.-',label="W from Dos")
# plt.plot(bins,Ave*W/bins,'.-',label="Ave from Dos * W/E")
# plt.plot(bins,Ave2*W/bins**2,'.-',label="Ave2 from Dos * $W^2/E^2$")
# plt.legend()
# plt.figure()

# plt.plot(bins,W,label="W")
# plt.plot(bins,calcW,label="calculated W")
# plt.legend()
# plt.figure()

plt.plot(bins,Ave,label="Ave")
plt.plot(bins,calcAve(bindex,A,B,C),label="calculated Ave")
plt.legend()
plt.figure()

plt.plot(bins,Ave2,label="Ave2")
plt.plot(bins,calcAve2(bindex,A,B,C),label="calculated Ave2")
plt.legend()

plt.show()




"""Extra"""


# calcW = np.zeros_like(bins)
# calcAve = np.copy(bins)
# calcAve2 = np.copy(bins)

# for x in range(0,len(bins)):
#     E0 = bins[x] - bstep*0.5
#     E1 = bins[x] + bstep*0.5
#     if C < 0:
#         F3 = B/(2*np.sqrt(-C))
#         F1 = special.erf(np.sqrt(-C)*E1-F3) - special.erf(np.sqrt(-C)*E0-F3)
#         calcW[x] = 0.5*np.sqrt(-np.pi/C)*np.exp(A+F3*F3)*F1
#         comp1 = np.exp(B*E0+C*E0*E0) - np.exp(B*E1+C*E1*E1)
#         comp2 = np.exp(B*E0+C*E0*E0)*(B-2*C*E0) - np.exp(B*E1+C*E1*E1)*(B-2*C*E1)
#         calcAve[x] = 1/(np.sqrt(-C*np.pi)*np.exp(F3**2))*comp1/F1+F3/np.sqrt(-C)
#         calcAve2[x] = 1/(2*(-C)**1.5*np.pi**0.5*np.exp(F3**2))*comp2/F1+(B**2-2*C)/(4*C**2)
#     else:
#         F3 = B/(2*np.sqrt(C))
#         F1 = special.erfi(np.sqrt(C)*E1+F3) - special.erfi(np.sqrt(C)*E0+F3)
#         calcW[x] = 0.5*np.sqrt(np.pi/C)*np.exp(A-F3*F3)*F1
#         comp1 = np.exp(B*E0+C*E0*E0) - np.exp(B*E1+C*E1*E1)
#         comp2 = np.exp(B*E0+C*E0*E0)*(B-2*C*E0) - np.exp(B*E1+C*E1*E1)*(B-2*C*E1)
#         calcAve[x] = -1/(np.sqrt(C*np.pi))*np.exp(F3**2)*comp1/F1-F3/np.sqrt(C)
#         calcAve2[x] = 1/(2*(C)**1.5*np.pi**0.5)*np.exp(F3**2)*comp2/F1+(B**2-2*C)/(4*C**2)
