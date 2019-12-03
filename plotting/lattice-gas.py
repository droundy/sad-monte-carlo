import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

with open('test.yaml','r') as stream:
    try:
        first_data_loaded = yaml.safe_load(stream)
    except yaml.yamlerror as exc:
        print(exc)

with open('partial_test.yaml','r') as stream:
    try:
        data_loaded = yaml.safe_load(stream)
    except yaml.yamlerror as exc:
        print(exc)
        
hist = data_loaded['bins']['histogram']

time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
energy_data = data_loaded['movies']['energy']
number_data = first_data_loaded['movies']['number']
#L = [1,2,3,4]
energy_resize = np.array(energy_data)
nlist = len(energy_data)
energy_resize.resize(8, 5)

fig, (ax0) = plt.subplots(1)
c = ax0.pcolor(energy_resize)
ax0.set_title('Lattice Gas Energy Pcolor')
fig.tight_layout()
plt.show()

#entropy_resize = np.array(entropy_data)
"""
flat_list = []
for sublist in entropy_data:
    print(sublist)
    for item in sublist:
        flat_list.append(item)
"""    
    
entropy_resize = []
for sublist in entropy_data[6:47]:
    #print(sublist)
    for item in sublist:
        #print(item)
        entropy_resize.append(item)
        
entropy_re = np.array(entropy_resize)
#print(entropy_re)
entropy_re.resize(47-6, 45)

#entropy_resize.resize(,)
#print(entropy_resize)

fig, (ax0) = plt.subplots(1)
c = ax0.pcolor(entropy_re)
ax0.set_title('Lattice Gas Entropy Pcolor')
fig.tight_layout()
plt.show()


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

#statistics


