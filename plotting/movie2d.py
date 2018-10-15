#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt
import glob

#~ def latex_float(x):
    #~ exp = np.log10(x*1.0)
    #~ if abs(exp) > 2:
        #~ x /= 10.0**exp
        #~ if ('%g' % x) == '1':
            #~ return r'10^{%.0f}' % (exp)
        #~ return r'%g\times 10^{%.0f}' % (x, exp)
    #~ else:
        #~ return '%g' % x

#~ allcolors = list(reversed(['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                           #~ 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']))

#~ my_energy = {}
#~ my_histogram = {}
#~ my_entropy = {}
#~ my_time = {}
#~ my_color = {}
#~ max_iter = 0
#~ my_gamma = {}
#~ my_gamma_t = {}
#~ Smin = None
#~ minT = 0.5
#~ fnames = sys.argv[1:]

for my_histogram in sorted(glob.iglob("samc-1e6-64-movie.yaml/h*.dat")):
    print(my_histogram)
