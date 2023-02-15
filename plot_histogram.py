import matplotlib.pyplot as plt
import numpy as np

def f(x,s,normed):
    if normed:
        return np.exp(-x**2/(2.0*s**2))/(s*np.sqrt(2.0*np.pi))
    else:
        return np.exp(-x**2/(2.0*s**2))

def getC(data):
    den, bin_edges = np.histogram(data[:,0], bins=1001, range=(-5.0,5.0))
    num, bin_edges = np.histogram(data[:,1], bins=1001, range=(-5.0,5.0))
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    return bin_centers, num/(den + 1e-10)
        

data = np.loadtxt('pairs.out')
data_V4a = np.loadtxt('pairs_V4a.out')

CF = getC(data)
CF_V4a = getC(data_V4a)

print(CF[0].shape,CF[1].shape)

plt.plot(CF[0], CF[1], 'r-')
plt.plot(CF_V4a[0], CF_V4a[1], 'b-')

xpts = np.linspace(-5,5,1001)
plt.plot(xpts, 1.0+2.75*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black')
plt.plot(xpts, 1.0+0.8*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black')
plt.plot(xpts, 1.0+0.4*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black')

plt.xlim([-0.25, 0.25])

#plt.hist(data[:,0], bins=1000, histtype='step', color='blue')
#plt.hist(data[:,1], bins=1000, histtype='step', color='red')

#xpts = np.linspace(-5,5,1001)
#plt.plot(xpts, 11500.0*f(xpts,np.sqrt(2.0), False), '-', color='black')
#plt.plot(xpts, 11500.0*f(xpts,1.0, False), '--', color='black')


plt.show()


# def f(x,s):
#     return np.exp(-x**2/(2.0*s**2))/(s*np.sqrt(2.0*np.pi))
# 
# #data = np.loadtxt('pairs_n100.out')
# #data2 = np.loadtxt('pairs_n1000.out')
# #data3 = np.loadtxt('pairs_n10000.out')
# 
# #data = np.loadtxt('pairs.out')
# #data = np.loadtxt('pairs_mult5.out')
# #data2 = np.loadtxt('pairs_mult5_nLoop100.out')
# data = np.loadtxt('pairs_mult3_nLoop100.out')
# data2 = np.loadtxt('pairs_mult3_nLoop100_FULL.out')
# data3 = np.loadtxt('pairs_mult3_nLoop1000_FULL.out')
# data3b = np.loadtxt('pairs_mult3_nLoop1000_FULL_nev5M.out')
# 
# print(data.shape)
# 
# plt.hist(data3b[:,0], bins=1000, density=True, histtype='step', color='blue')
# plt.hist(data[:,1], bins=1000, density=True, histtype='step', color='red')
# plt.hist(data2[:,1], bins=1000, density=True, histtype='step', color='green')
# plt.hist(data3[:,1], bins=1000, density=True, histtype='step', color='purple')
# plt.hist(data3b[:,1], bins=1000, density=True, histtype='step', color='cyan')
# 
# xpts = np.linspace(-5,5,1001)
# plt.plot(xpts, f(xpts,1.0), '-', color='black')
# plt.plot(xpts, f(xpts,1.0)*(1.0+np.exp(-0.5*xpts**2*(5.0/0.19733)**2)), '--', color='black')
# plt.plot(xpts, f(xpts,1.0)*(1.0+0.95*np.exp(-0.5*xpts**2*(5.0/0.19733)**2)), ':', color='black')
# plt.plot(xpts, f(xpts,np.sqrt(2.0)), '-', color='black')

#plt.show()
#plt.savefig('histogram.png')
