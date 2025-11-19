import numpy as np
import matplotlib.pyplot as plt

deltanll = np.loadtxt('deltanll.txt')
nll = np.loadtxt('nll.txt')
likelihood = np.loadtxt('likelihood.txt')


plt.plot(likelihood[:,0],likelihood[:,1])
plt.savefig('likelihood.png')
plt.clf()

plt.plot(nll[:,0],nll[:,1])
plt.savefig('nll.png')
plt.clf()

plt.plot(deltanll[:,0],deltanll[:,1])
plt.savefig('deltanll.png')
plt.clf()