import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(0, 1, 10000)

a = 0.1
b = 2*a

e1 = 1
e2 = 0.5

E = np.zeros_like(r)

E = -e2*np.exp(-(r-b)**2/(b-a)**2)
E[r<b] = (r[r<b]-b)**2/((b-a)**2)*e2 - e2

E = (r-b)**2/((b-a)**2)*e2 - e2

E[r<a] = r[r<a]**2/a**2 - e1

plt.plot(r, E)

plt.show()