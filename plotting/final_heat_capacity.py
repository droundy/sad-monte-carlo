import sys, re
import numpy as np
import matplotlib.pyplot as plt

# HOW TO RUN:
# python3 plotting/final_heat_capacity.py lj-sad-31-minT001-de001.yaml lj-sad*yaml lj-samc-1e*yaml 0.01

def latex_float(x):
    exp = int(np.log10(x*1.0))
    if abs(exp) > 2:
        x /= 10.0**exp
        if ('%.1g' % x) == '1':
            return r'10^{%.0f}' % (exp)
        return r'%.1g\times10^{%.0f}' % (x, exp)
    else:
        return '%g' % x

def lookup_entry(entry, yaml_data):
    x = re.search('{}: (\S+)'.format(entry), yaml_data)
    if x:
        return float(x.group(1))

def fix_fname(fname):
    if fname[-5:] == '.yaml':
        return fname[:-5]
    return fname

# Import the data for
ref1 = 'LJ31_Cv_Reference'
data1 = np.loadtxt(ref1 + '.csv', delimiter = ',', unpack = True)
ref1_T = data1[0]
ref1_heat_capacity = (data1[1]-3/2)*31

ref2 = 'LJ31_Cv_Reference_alt'
data2 = np.loadtxt(ref2 + '.csv', delimiter = ',', unpack = True)
ref2_T = data2[0]
ref2_heat_capacity = data2[1]

ref3 = 'LJ31_Cv_Reference_3'
data3 = np.loadtxt(ref3 + '.csv', delimiter = ',', unpack = True)
ref3_T = data3[0]
ref3_heat_capacity = (data3[1]-3/2)*31

ref4 = 'LJ31_Cv_Reference_4'
data4 = np.loadtxt(ref4 + '.csv', delimiter = ',', unpack = True)
ref4_T = data4[0]
ref4_heat_capacity = (data4[1]-3/2)*31

my_energy = {}
my_de = {}
my_histogram = {}
my_entropy = {}
my_time = {}
my_color = {}
max_iter = 0
Smin = None
my_minT = {}
my_too_lo = {}
my_too_hi = {}
minT = None

allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
                           'xkcd:lightblue', 'xkcd:puke', 'xkcd:puce', 'xkcd:turquoise']))

# will be None if a float i.e. the user input a min temperature.
last_fname = re.search('[a-zA-Z]', sys.argv[-1])

if not last_fname:
    fnames = [fix_fname(f) for f in sys.argv[1:-1]]
else:
    fnames = [fix_fname(f) for f in sys.argv[1:]]


def heat_capacity(T, E, S):
    C = np.zeros_like(T)
    for i in range(len(T)):
        boltz_arg = S - E/T[i]
        P = np.exp(boltz_arg - boltz_arg.max())
        P = P/P.sum()
        U = (E*P).sum()
        C[i] = ((E-U)**2*P).sum()/T[i]**2
    return C

fig, ax = plt.subplots(figsize=[5, 4])

# inset axes....
axins = ax.inset_axes([0.15, 0.50, 0.40, 0.45])

for fname in fnames:
    print(fname)
    if len(sys.argv) > 1 and not last_fname: # the user has input a min temp.
        minT = float(sys.argv[-1])
        print('minT =', minT)
    else:
        for i in range(len(fname)-5):
            if fname[i:i+len('minT')] == 'minT':
                minT = float(fname[i+len('minT'):].split('-')[0])
                print('minT =', minT)

    my_energy[fname] = np.loadtxt(fname+'.energy')
    my_entropy[fname] = np.loadtxt(fname+'.entropy')
    my_color[fname] = allcolors.pop()

    T = np.linspace(minT,0.4,500)
    ax.plot(T, heat_capacity(T, my_energy[fname], my_entropy[fname][-1:]), my_color[fname],
        label=fname, linewidth=2)
    axins.plot(T, heat_capacity(T, my_energy[fname], my_entropy[fname][-1:]), my_color[fname],
        label=fname, linewidth=2)

plt.ylabel('heat capacity')
plt.xlabel('temperature')

# REM (Replica exchange Monte Carlo) REFERENCE comes from https://pubs.acs.org/doi/full/10.1021/jp072929c
# They use a radius of Rc = ? and an energy bin size ?

# RESTMC REFERENCE comes from Replica exchange statistical temperature Monte Carlo by Kim, Keyes, and Straub
# They use a radius of Rc = 2.5/sigma and an energy bin size 0.2 at minimum

# MCET (Monte Carlo Entropic Tempering) REFERENCE comes from https://journals.aps.org/pre/pdf/10.1103/PhysRevE.63.010902
# They use a radius of Rc = ? and an energy bin size ?

# EMC (Exchange Monte Carlo) REFERENCE comes from https://aip.scitation.org/doi/pdf/10.1063/1.2202312?class=pdf&
# The referee directed us to this paper!

plt.plot(ref1_T, ref1_heat_capacity, '-', color='black', label='REM REFERENCE',linewidth=1)
plt.plot(ref2_T, ref2_heat_capacity, '-', color='grey', label='RESTMC REFERENCE',linewidth=1)
plt.plot(ref3_T, ref3_heat_capacity, '--', color='blue', label='MCET REFERENCE',linewidth=1)
plt.plot(ref4_T, ref4_heat_capacity, '--', color='black', label='ECM REFERENCE',linewidth=2)

axins.plot(ref1_T, ref1_heat_capacity, '-', color='black', label='REM REFERENCE',linewidth=1)
axins.plot(ref2_T, ref2_heat_capacity, '-', color='grey', label='RESTMC REFERENCE',linewidth=1)
axins.plot(ref3_T, ref3_heat_capacity, '--', color='blue', label='MCET REFERENCE',linewidth=1)
axins.plot(ref4_T, ref4_heat_capacity, '--', color='black', label='ECM REFERENCE',linewidth=2)

ax.legend(loc='best')
plt.xlim(0.01,0.4)
plt.ylim(40,140)

# sub region of the original image
x1, x2, y1, y2 = 0.01, 0.04, 40, 100
dx, dy = 0.005, 10

axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels(np.around(np.arange(x1,x2+dx,dx),3).tolist())
axins.set_yticklabels(np.around(np.arange(y1,y2+dy,dy),3).tolist())
#axins.set_major_formatter(mtick.FormatStrFormatter('%.3'))
ax.indicate_inset_zoom(axins)


plt.show()
