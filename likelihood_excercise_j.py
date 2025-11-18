import numpy as np
import scipy
import matplotlib.pyplot as plt

def gaussian(x, mean, sigma):
    return np.exp(-(x - mean)**2 / (2 * sigma**2))/(np.sqrt(2*np.pi) * sigma)

def poisson(k, lam):
    """Poisson PMF model (not normalized to data size)."""
    return scipy.stats.poisson.pmf(k, lam)

Nsr = 250
Ncr = 500
B1mc = 100
b1mc_sig = 2
Tau2mc = 5
tau2mc_sig = 0.2


def likelihood(x,S,Nsr,Ncr,B1mc,b1mc_sig,Tau2mc,tau2mc_sig):
    l = poisson(int(round(S+x[0]+x[1])),Nsr)*poisson(int(round(x[1]*x[2])),Ncr)*gaussian(x[2],Tau2mc,tau2mc_sig)*gaussian(x[0],B1mc,b1mc_sig)
    return -2*np.log(l)



solutions = np.empty((0,3))
fun = np.array([])
for S in range(100):
    x0 = [100,100,5]
    res = scipy.optimize.minimize(likelihood, x0,args=(S,Nsr,Ncr,B1mc,b1mc_sig,Tau2mc,tau2mc_sig), method="Nelder-Mead")
    solutions = np.append(solutions,res.x)
    fun = np.append(fun,res.fun)

plt.plot(np.arange(100),fun)
plt.show()
