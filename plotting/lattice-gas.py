import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

with open('test.yaml','r') as stream:
    try:
        data_loaded = yaml.safe_load(stream)
    except yaml.yamlerror as exc:
        print(exc)
        
hist = data_loaded['bins']['histogram']

time_frame = data_loaded['movies']['time']
print('entropy_data', data_loaded['movies']['entropy'])
print('entropy_data lens', [len(x) for x in data_loaded['movies']['entropy']])
entropy_data = data_loaded['movies']['entropy'][7:]
moves = data_loaded['movies']['time'][7:]

energy_data = data_loaded['movies']['energy']
number_data = np.array(data_loaded['movies']['number'])

energy_resize = np.array(energy_data)
print('energy_resize size', energy_resize.shape)
nlist = len(energy_data)
energy_resize.resize(8, 5)
number_data.resize(8, 5)

fig, (ax0) = plt.subplots(1)
c = ax0.pcolor(energy_resize)
ax0.set_title('Lattice Gas Energy Pcolor')
fig.tight_layout()

fig, (ax0) = plt.subplots(1)
c = ax0.pcolor(number_data)
ax0.set_title('Lattice Gas Number Pcolor')
fig.tight_layout()

"""
flat_list = []
for sublist in entropy_data:
    print(sublist)
    for item in sublist:
        flat_list.append(item)
"""    

plt.figure()
for t in range(len(entropy_data)):
    print('time', moves[t])
    S = np.array(entropy_data[t])
    S.resize(8,5)
    plt.title(f'{moves[t]} moves')
    plt.pcolor(number_data, energy_resize, S)
    plt.pause(1)

entropy_resize = []
for sublist in entropy_data[6:47]:
    for item in sublist:
        entropy_resize.append(item)

entropy_re = np.array(entropy_resize)
entropy_re.resize(47-6, 45)

#matplotlib.pyplot.pcolor(*args, alpha=None, norm=None, cmap=None, vmin=None, vmax=None, data=None, **kwargs)

#entropy_resize.resize(,)
#print(entropy_resize)
fig, (ax0) = plt.subplots(1)
c = ax0.pcolor(entropy_re)
ax0.set_title('Lattice Gas Entropy Pcolor')
fig.tight_layout()


#graphing
"""
plt.plot(hist)
plt.ylabel('Histogram of Entropy')

#plt.plot(time_frame, energy_data)
#plt.label('Histogram of Entropy')

plt.plot(time_frame[0:45], energy_data[0:45])
plt.title('Energy vs Time of Lattice Gas')
plt.ylabel('energy')
plt.xlabel('time')
#plt.legend(loc = 'best')
#plt.plot(time_frame, entropy_data[0:44])
"""
plt.show()

