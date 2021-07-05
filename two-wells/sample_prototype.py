import numpy as np
import scipy.special as spl
import scipy.interpolate as pol
import matplotlib.pyplot as plt
import math as mth


def V(n):
    return np.pi**(n/2)/(spl.gamma(n/2+1))


def pdf(x, dim):
    return  (1-x**2)**(dim/2) * V(dim)/V(dim+1)


def cdf(x, dim, du):
    val = 0
    for u in np.arange(-1+du,x-du,du):
        val += (pdf(u, dim) + 4 * pdf((2*u + du) / 2, dim) + pdf(u+du, dim)) * du/6
    return val




def generate_stencil(num_points, dim):
    numerical_precision_mult = 100
    sten = np.zeros(num_points*(dim-1))
    us = np.linspace(-1,1,(num_points-1) * numerical_precision_mult)
    du = us[1]-us[0]
    
    for N in range(0, dim-1):
        val = 0
        w = 1
        for i in range(1,(num_points-1)*numerical_precision_mult):
            if(i % numerical_precision_mult == 0):
                sten[N*num_points + w] = val
                w+=1

            val += (pdf(us[i-1], dim-N-1) + 4 * pdf((us[i] + us[i-1]) / 2, dim-N-1) + pdf(us[i], dim-N-1)) * du/6

        sten[N*num_points + w] = 1.#val

    return sten


def find_bin(data, lower_index, upper_index, val):
    if(data[upper_index] <= val):
        return [data[upper_index-1], data[upper_index], upper_index - lower_index]

    first_lower = lower_index
    guess = (lower_index + upper_index)//2

    while(lower_index != upper_index - 1):
        if(data[guess] <= val):
            lower_index = guess
        elif(data[guess] > val):
            upper_index = guess
        else:
            print("Something's wrong in the search?")

        guess = (lower_index + upper_index)//2

    return [data[lower_index], data[upper_index], lower_index - first_lower]


class BallInvCdf:
    def __init__(self, num_points, dim):
        self.num_points = num_points
        self.dim = dim
        self.dx = 2/(num_points)

        self.stencils = np.zeros(self.num_points*(self.dim-1))
        numerical_precision_mult = 100

        us = np.linspace(-1,1,(num_points-1) * numerical_precision_mult)
        du = us[1]-us[0]
    
        for N in range(0, self.dim-1):
            val = 0
            w = 1
            for i in range(1,(self.num_points-1)*numerical_precision_mult):
                if(i % numerical_precision_mult == 0):
                    self.stencils[N*self.num_points + w] = val
                    w+=1

                val += (pdf(us[i-1], self.dim-N-1) + 4 * pdf((us[i] + us[i-1]) / 2, self.dim-N-1) + pdf(us[i], self.dim-N-1)) * du/6

            self.stencils[N*self.num_points + w] = 1.#val



    def eval(self, probability, which_dim):
        L, U, j = find_bin(self.stencils, which_dim*self.num_points, (which_dim+1)*self.num_points-1, probability)
        slope = self.dx/(U-L)
        return slope * (probability - L) + -1 + j*self.dx

    def print_data(self):
        print(self.stencils)

    def sample(self, radius):
        R = radius
        sample = np.zeros(self.dim)
        for i in range(0,self.dim-1):
            sample[i] = self.eval(np.random.uniform(0,1), i)
            sample[i] *= R
            R = np.sqrt(R**2 - sample[i]**2)

        sample[self.dim - 1] = R * np.random.uniform(-1,1)

        return sample




class SystemInvCdf:
    def __init__(self, num_points, dim, r1, r2):
        self.num_points = num_points
        self.dim = dim
        self.dx = 2/(num_points)
        self.r1 = r1
        self.r2 = r2
        x1 = np.linspace(-r1, r1+2*r2, self.num_points)
        self.dx1_ball1 = x1[1] - x1[0]

        ######## The first axis, ball 1 ########

        self.stencils = np.zeros(self.num_points*(self.dim))
        numerical_precision_mult = 100

        val = 0
        for w in range(self.num_points-1):
            us = np.linspace(x1[w], x1[w+1], numerical_precision_mult)
            du = us[1] - us[0]
            for i in range(numerical_precision_mult-1):
                # use midpoint method
                val += du*self.pdf_x1_nonnormalized(0.5*(us[i+1]+us[i]))
            self.stencils[w+1] = val

        self.V_ball1 = val

        for i in range(0, self.num_points):
            self.stencils[i] /= val

        ######## Find the volume of the cylinder, second ball, and total volume ########

        self.V_ball2 = V(dim) * (r2)**dim / 2.
        self.V_cyl = (r1 + 2*r2 - np.sqrt(r1**2 - r2**2)) * V(dim-1) * (r2)**dim
        self.V_tot = self.V_ball1 + self.V_ball2 + self.V_cyl

        ######## calculate Probabilities of being in Ball1, the cylinder, or Ball 2 ########

        self.P_ball1 = self.V_ball1 / self.V_tot
        self.P_cyl = self.V_cyl / self.V_tot
        self.P_ball2 = self.V_ball2 / self.V_tot
        
        ######## The rest follow a uniform sphere ########

        us = np.linspace(-1,1,(num_points-1) * numerical_precision_mult)
        du = us[1]-us[0]
    
        for N in range(1, self.dim):
            val = 0
            w = 1
            for i in range(1,(self.num_points-1)*numerical_precision_mult):
                if(i % numerical_precision_mult == 0):
                    self.stencils[N*self.num_points + w] = val
                    w+=1

                val += (pdf(us[i-1], self.dim-N) + 4 * pdf((us[i] + us[i-1]) / 2, self.dim-N) + pdf(us[i], self.dim-N)) * du/6

            self.stencils[N*self.num_points + w] = 1.#val



    def pdf_x1_nonnormalized(self, x1):#NOT NORMALIZED
        if(x1 <= np.sqrt(self.r1**2 - self.r2**2)):
            return (self.r1**2-x1**2)**((self.dim-1)/2) * V(self.dim - 1)

        elif(x1 < self.r1+self.r2):
            return (self.r2)**(self.dim-1) * V(self.dim - 1)

        else:
            return (self.r2**2-(x1-self.r1-self.r2)**2)**((self.dim-1)/2) * V(self.dim - 1)


    def pdf_x1(self, x1):
        if(x1 <= np.sqrt(self.r1**2 - self.r2**2)):
            return (self.r1**2-x1**2)**((self.dim-1)/2) * V(self.dim - 1) / self.V_tot

        elif(x1 < self.r1+self.r2):
            return (self.r2)**(self.dim-1) * V(self.dim - 1) / self.V_tot

        else:
            return (self.r2**2-(x1-self.r1-self.r2)**2)**((self.dim-1)/2) * V(self.dim - 1) / self.V_tot




    def eval(self, probability, which_dim):
        L, U, j = find_bin(self.stencils, which_dim*self.num_points, (which_dim+1)*self.num_points-1, probability)
        if(which_dim == 0):
            slope = self.dx1_ball1/(U-L)
            return slope * (probability - L) + -self.r1 + j*self.dx1_ball1
        # elif(which_dim == 1):
        #     slope = self.dx1_ball2/(U-L)
        #     return slope * (probability - L) + self.r1 + self.r2 + j*self.dx1_ball2
        else:
            slope = self.dx/(U-L)
        return slope * (probability - L) + -1. + j*self.dx


    def print_data(self, which_dim):
        print(self.stencils[which_dim * self.num_points : (which_dim+1) * self.num_points])


    def sample(self):
        sample = np.zeros(self.dim)
        ### Determine if in ball 1 ###
        random_num = np.random.uniform(0,1)

        x1 = self.eval(np.random.uniform(0,1), 0)
        sample[0] = x1
        if x1 <= np.sqrt(self.r1**2 - self.r2**2):
            # we are in the first ball
            R = np.sqrt(self.r1**2 - x1**2)
        elif x1 <self.r1+self.r2:
            # we're in the cylinder
            R = self.r2
        else:
            # We're in the final hemisphere
            R = np.sqrt(self.r2**2 - (x1-self.r1-self.r2)**2)

        for i in range(2,self.dim):
            sample[i-1] = self.eval(np.random.uniform(0,1), i)
            sample[i-1] *= R
            R = np.sqrt(R**2 - sample[i-1]**2)

        sample[self.dim - 1] = R * np.random.uniform(-1,1)

        return sample

    


# def lin_interp(data, lower_index, upper_index, dx, eval):
#     L, U, j = find_bin(data, lower_index, upper_index, eval)
#     slope = dx/(U-L)
#     return slope * (eval - L) + -1 + j*dx





# dx = 0.00001
# a = []
# for i in np.arange(0,1,dx):
#     a.append(i)

# y = []
# ex = np.arange(0,1,0.1)
# for i in range(0,len(ex)):
#     y.append(lin_interp(a, 0, len(a)-1, dx, ex[i]**2))

# plt.plot(y,ex)


# f = BallInvCdf(100,2)

# f.print_data()



# ex = np.linspace(0.,1.,200)
# ys = []

# for x in ex:
#     ys.append(f.eval(x, 0))

# plt.plot(ex,ys)


# print(f.eval(1.,0))

# plt.figure()

# t = np.linspace(0,2*np.pi, 100)
# plt.plot(2*np.cos(t), 2*np.sin(t),'--')

# for i in range(0,500):
#     s = f.sample(2)
#     plt.plot(s[0],s[1],'.')


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')


# f = BallInvCdf(100,3)

# for i in range(0,500):
#     s = f.sample(1.5)
#     ax.scatter(s[0], s[1], s[2])


# plt.figure()

g = SystemInvCdf(10, 3, 0.75, 0.5)
g.print_data(1)
print("\n\nVolume of Ball 1: " + str(g.V_ball1))
print("Volume of Cylinder: " + str(g.V_cyl))
print("Volume of Ball 2: " + str(g.V_ball2))
print("Total Volume: " + str(g.V_tot))
print("Probability of being in Ball 1: " + str(g.P_ball1))
print("Probability of being in Cylinder: " + str(g.P_cyl))
print("Probability of being in Ball 2: " + str(g.P_ball2))


ex = np.linspace(-0.75, 0.75 + 2*0.5, 1000)
p = np.zeros_like(ex)
for i in range(len(p)):
    p[i] = g.pdf_x1(ex[i])

hist = np.zeros_like(p)
N_hist = 1000000
for percent in range(100):
    for _ in range(N_hist//100):
        s = g.sample()
        hist[abs(ex-s[0]).argmin()] += 1
    print(f'{percent}% done')
hist *= p.sum()/hist.sum() # normalize

plt.plot(ex,p)
plt.plot(ex, hist, '.')
plt.title('pdf inside the big well')

plt.figure('x1')


ex = np.linspace(0, 1, 1000)
p = []
for x in ex:
    p.append(g.eval(x, 0))

plt.plot(ex,p)
plt.title('inv cdf in the big well')

# plt.figure()

# num_samples = 5000

# samples = np.zeros((2,num_samples))

# for i in range(0,num_samples):
#     samples[:,i] = g.sample()

# plt.plot(samples[0,:], samples[1,:],'.')

g = SystemInvCdf(1024, 2, 0.75, 0.5)

fig = plt.figure()
# ax = fig.add_subplot(projection='3d')


num_samples = 10000

samples = np.zeros((2,num_samples))

for i in range(0,num_samples):
    samples[:,i] = g.sample()

plt.plot(samples[0,:], samples[1,:],'.')
# ax.scatter(samples[0,:], samples[1,:], samples[2,:])


plt.show()