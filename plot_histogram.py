import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sf

nBins = 1001

def f(x,s,normed):
    if normed:
        return np.exp(-x**2/(2.0*s**2))/(s*np.sqrt(2.0*np.pi))
    else:
        return np.exp(-x**2/(2.0*s**2))
    
def getavg(a, b, x0, x1):
    return np.mean(b[((a>=x0) & (a<=x1)) | ((a<=-x0) & (a>=-x1))])

def geterr(x, y, nev, npairs):
    x.reshape(nev,npairs)

def getC(data, symmetric=True):
    npairs = 10
    nev = len(data[:,0])//npairs
    den, bin_edges = np.histogram(data[:,0], bins=nBins, range=(-5.0,5.0))
    num, bin_edges = np.histogram(data[:,1], bins=nBins, range=(-5.0,5.0))
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    if symmetric:
        return bin_centers, 0.5*(num+num[::-1])/(0.5*(den+den[::-1]) + 1e-10)
    else:
        return bin_centers, num/(den + 1e-10)
    

def getden(data, symmetric=True):
    den, bin_edges = np.histogram(data[:,0], bins=nBins, range=(-5.0,5.0))
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    if symmetric:
        return bin_centers, 0.5*(den+den[::-1])
    else:
        return bin_centers, den
    
def getnum(data, symmetric=True):
    num, bin_edges = np.histogram(data[:,1], bins=nBins, range=(-5.0,5.0))
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    if symmetric:
        return bin_centers, 0.5*(num+num[::-1])
    else:
        return bin_centers, num
    

def F(x):
    a = 5.0/(np.sqrt(2.0)*0.19733)
    return 0.5*(1.0+sf.erf(a*x))

def KS(x):
    return 3.0-2.0*np.maximum(F(x), 1.0-F(x))
    
        
plotCF = True
logPlot = False
differential = False

data = np.loadtxt('n4_TRIAL4b_nL100.out')

#data_all = np.loadtxt('pairs_n3_n15000000_allorder.out')
#data2_all = np.loadtxt('pairs_n3_n15000000_allorder_run2.out')
#data3_all = np.loadtxt('pairs_n3_n15000000_allorder_run3.out')
#data_2nd = np.loadtxt('pairs_n3_n15000000_2ndorder.out')
#data2_2nd = np.loadtxt('pairs_n3_n15000000_2ndorder_run2.out')
#data3_2nd = np.loadtxt('pairs_n3_n15000000_2ndorder_run3.out')

#data = np.concatenate((data_all,data2_all,data3_all))
#data2 = np.concatenate((data_2nd,data2_2nd,data3_2nd))


if plotCF:
    CF = getC(data, False)
    #CF2 = getC(data2, False)
    ##CF3 = getC(data3, False)
    ##CF4 = getC(data4, False)

    print(CF[0].shape,CF[1].shape)

    if logPlot:
        CF0 = (CF[0])[CF[0]>0], (CF[1])[CF[0]>0]
        plt.semilogy(CF0[0]**2, CF0[1]-1.0, 'r-')
        #plt.plot(CF_V4a[0], CF_V4a[1], 'b-')

        x2pts = np.linspace(0,5**2,nBins)
        plt.semilogy(x2pts, np.exp(-0.5*x2pts*(5.0/0.19733)**2), '-', color='black')
        #plt.semilogy(x2pts, 2.36031*np.exp(-0.5*x2pts*(5.0/0.19733)**2), '--', color='black')
        plt.plot(x2pts, 1.41*np.exp(-0.5*x2pts*(5.0/0.19733)**2), ':', color='black')

        plt.xlim([0.0, 0.125**2])
        plt.ylim([3e-2, 2.1])
    else:
        if differential:
            print(getavg(CF[0], CF[1], 0.15, 0.5))
            plt.plot(CF[0], CF[1]/getavg(CF[0], CF[1], 0.2, 0.5)-2.36031*np.exp(-0.5*CF[0]**2*(5.0/0.19733)**2), 'ro-')

            xpts = np.linspace(-5,5,nBins)
            plt.plot(xpts, 0.0*xpts+1.0, ':', color='black')

            plt.xlim([-0.5, 0.5])
            plt.ylim([0.9, 2.5])
        else:
            print(getavg(CF[0], CF[1], 0.2, 0.5))
            #print(getavg(CF2[0], CF2[1], 0.15, 0.5))
            #print(getavg(CF3[0], CF3[1], 0.15, 0.5))
            #print(getavg(CF4[0], CF4[1], 0.15, 0.5))
            plt.plot(CF[0], CF[1]/getavg(CF[0], CF[1], 0.2, 0.5), 'r-')
            #plt.plot(CF2[0], CF2[1]/getavg(CF2[0], CF2[1], 0.2, 0.5), 'b-')
            #plt.plot(CF3[0], CF3[1]/getavg(CF3[0], CF3[1], 0.2, 0.5), 'go-')
            #plt.plot(CF4[0], CF4[1]/getavg(CF4[0], CF4[1], 0.2, 0.5), 'o-', color='purple')

            xpts = np.linspace(-5,5,nBins)
            plt.plot(xpts, 1.0+np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black') #n=2
            #plt.plot(xpts, 1.0+1.0384*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '--', color='black') #n=3
            #plt.plot(xpts, 1.0+1.0322*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), ':', color='black') #n=3
            #plt.plot(xpts, 1.0+1.09567*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), ':', color='black') #n=4
            #plt.plot(xpts, 1.0+1.17695*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black') #n=5
            #plt.plot(xpts, 1.0+2.36031*np.exp(-0.5*xpts**2*(5.0/0.19733)**2), '-', color='black') #n=10
            plt.plot(xpts, 1.0+(2.0/3.0)*np.exp(-0.5*xpts**2*(5/0.19733)**2), '--', color='black') #n=10
            plt.plot(xpts, 1.0+(2.0/4.0)*np.exp(-0.5*xpts**2*(5/0.19733)**2), ':', color='black') #n=10
            plt.plot(xpts, 1.0+1.15*np.exp(-0.5*xpts**2*(5/0.19733)**2), ':', color='black') #n=10
            #plt.plot(xpts, KS(xpts), '.-', color='black') #n=10
            plt.plot(xpts, np.exp(-0.5*xpts**2), ':', color='black') #n=10
            plt.plot(xpts, np.exp(-0.5*xpts**2*0.25**2), ':', color='black') #n=10


            plt.xlim([-0.2, 0.2])
            plt.ylim([0.95, 2.05])

else:
    num = getnum(data, False)
    den = getden(data, False)
    num2 = getnum(data2, False)
    den2 = getden(data2, False)

    #plt.plot(num[0], num[1], 'r-')
    #plt.plot(den[0], den[1], 'b-')
    plt.plot(num2[0], num2[1], 'g-')
    plt.plot(den2[0], den2[1], '-', color='purple')

    xpts = np.linspace(-5,5,nBins)
    norm, width, norm2 = 110, 0.15, 615
    #plt.plot(xpts, norm*np.exp(-width*xpts**2), '--', color='black')
    #plt.plot(xpts, norm*np.exp(-width*xpts**2)*(1.0+np.exp(-0.5*xpts**2*(5/0.19733)**2)), '-', color='black')
    plt.plot(xpts, norm2*(1.0+np.exp(-0.5*xpts**2*(5/0.19733)**2)), '-', color='black')
    plt.plot(xpts, norm2*(1.0+0.75*np.exp(-0.5*xpts**2*(5/0.19733)**2)), '--', color='black')

    plt.xlim([-5.0, 5.0])
    #plt.ylim([0.0, 2.5])



#plt.hist(data[:,0], bins=nBins, histtype='step', color='blue')
#plt.hist(data[:,1], bins=nBins, histtype='step', color='red')

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
