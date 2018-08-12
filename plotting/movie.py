#!/usr/bin/python3

import yaml, sys
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

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
    energy = np.array(data['movies']['energy'])
    time = np.array(data['movies']['time'])
    S = np.array(data['movies']['entropy'])
    histogram = np.array(data['movies']['histogram'])

    print(energy.shape, S.shape)

    for i,t in enumerate(time):
        plt.figure('Entropy')
        plt.plot(energy, S[i,:], label='t=%g'%t)
        plt.title('$t=%.3g/%.3g' % (t, time[-1]))
        plt.ylabel('$S$')
        plt.legend(loc='best')

        plt.figure('Histogram')
        plt.title('$t=%.3g/%.3g' % (t, time[-1]))
        plt.ylabel('histogram')
        plt.plot(energy, histogram[i,:], label='t=%g'%t)
        #plt.legend(loc='best')
        plt.pause(5.0)
        print('on time', t)

plt.ioff()
plt.show()
