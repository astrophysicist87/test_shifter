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
        

logPlot = False

data = np.loadtxt('pairs.out')
#data2 = np.loadtxt('pairs_n10_wNorm_nL100.out')
#data3 = np.loadtxt('pairs_n3_nL100.out')
#data4 = np.loadtxt('pairs_n3FULL_nL100.out')

CF = getC(data)
#CF2 = getC(data2)
#CF3 = getC(data3)
#CF4 = getC(data4)

print(CF[0].shape,CF[1].shape)

if logPlot:
    CF0 = (CF[0])[CF[0]>0], (CF[1])[CF[0]>0]
    plt.semilogy(CF0[0]**2, CF0[1]-1.0, 'r-')
    #plt.plot(CF_V4a[0], CF_V4a[1], 'b-')

    x2pts = np.linspace(0,5**2,1001)
    plt.semilogy(x2pts, np.exp(-0.5*x2pts*(5.0/0.19733)**2), '-', color='black')
    plt.semilogy(x2pts, 0.91*np.exp(-0.5*x2pts*(5.0/0.19733)**2), '--', color='black')
    #plt.plot(xpts, 1.0+0.4*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), ':', color='black')

    plt.xlim([0.0, 0.09**2])
    plt.ylim([5e-2, 1.1])
else:
    plt.plot(CF[0], CF[1], 'r-')
    #plt.plot(CF2[0], CF2[1], 'b-')
    #plt.plot(CF3[0], CF3[1], 'g-')
    #plt.plot(CF4[0], CF4[1], '-', color='purple')

    xpts = np.linspace(-5,5,1001)
    plt.plot(xpts, 1.0+np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black')
    #plt.plot(xpts, 1.0+0.91*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '--', color='black')
    #plt.plot(xpts, 1.0+0.4*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), ':', color='black')

    plt.xlim([-0.25, 0.25])
    plt.ylim([0.9, 2.2])

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
