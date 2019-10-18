import sys, re
import numpy as np
import matplotlib.pyplot as plt

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
dirname = 'LJ31_Cv_Reference'
data = np.loadtxt(dirname + '.csv', delimiter = ',', unpack = True)
ref_T = data[0]
ref_heat_capacity = (data[1]-3/2)*31
print(ref_heat_capacity)


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
    plt.plot(T, heat_capacity(T, my_energy[fname], my_entropy[fname][-1:]), my_color[fname],
        label=fname)

plt.ylabel('heat capacity')
plt.xlabel('temperature')
plt.plot(ref_T, ref_heat_capacity, '-', color='black', label='Reference Cv',linewidth=1)
plt.legend(loc='best')
plt.show()




