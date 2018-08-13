#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt

allcolors = ['g','r','b',]

my_energy = {}
my_histogram = {}
my_entropy = {}
my_time = {}
my_color = {}
max_iter = 0
for fname in sys.argv[1:]:
    print(fname)
    with open(fname) as f:
        yaml_data = f.read()
    data = yaml.load(yaml_data)
    data['bins']['histogram'] = np.array(data['bins']['histogram'])
    data['bins']['lnw'] = np.array(data['bins']['lnw'])
    print('bins:', data['bins'].keys())
    print('movies:', data['movies'].keys())
    print(data.keys())
    my_color[fname] = allcolors.pop()
    my_energy[fname] = np.array(data['movies']['energy'])
    my_time[fname] = np.array(data['movies']['time'])
    if len(my_time[fname]) > max_iter:
        max_iter = len(my_time[fname])
    my_entropy[fname] = np.array(data['movies']['entropy'])
    my_histogram[fname] = np.array(data['movies']['histogram'])

plt.ion()
while True:
    for i in range(max_iter):
        plt.figure('Histogram')
        plt.cla()
        plt.figure('Entropy')
        plt.cla()
        for fname in my_energy.keys():
            if i < len(my_time[fname]):
                t = my_time[fname][i]
                plt.figure('Entropy')
                if i == 0:
                    plt.plot(my_energy[fname], my_entropy[fname][i,:], my_color[fname],
                             label=fname)
                else:
                    plt.plot(my_energy[fname], my_entropy[fname][i-1,:], my_color[fname],
                             alpha=0.5)
                    plt.plot(my_energy[fname], my_entropy[fname][i,:], my_color[fname])
                plt.title('$t=%.3g/%.3g' % (t, my_time[fname][-1]))
                plt.ylabel('$S$')
                plt.legend(loc='best')
                plt.figure('Histogram')
                plt.title('$t=%.3g/%.3g' % (t, my_time[fname][-1]))
                plt.ylabel('histogram')
                if i == 0:
                    plt.plot(my_energy[fname], my_histogram[fname][i,:], my_color[fname],
                             label=fname)
                else:
                    plt.plot(my_energy[fname], my_histogram[fname][i-1,:], my_color[fname],
                             alpha=0.5)
                    plt.plot(my_energy[fname], my_histogram[fname][i,:], my_color[fname])
                plt.legend(loc='best')
        plt.pause(1.0)
        print('on time', t)

plt.ioff()
plt.show()
