import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker
T = np.linspace(320, 900, num=50)
c_surf = np.logspace(20, 22, num=50, endpoint=True)

X, Y = np.meshgrid(T, c_surf)

nu_i = 1e13*np.exp(-1.0/8.617e-5/X)
nu_m = 1e-7*np.exp(-0.2/8.617e-5/X)/(6.0*6.3e28*1.1e-10*1.1e-10)

R_trap = 1/(1+nu_i/nu_m/Y)

levels = np.logspace(
    np.log10(np.min(R_trap)),
    np.log10(np.max(R_trap)),
    1000)
levels2 = np.logspace(
    np.log10(np.min(R_trap)),
    np.log10(np.max(R_trap)),
    10)
locator = ticker.LogLocator(base=10)
plt.yscale("log")
CS = plt.contourf(X, Y, R_trap, locator=locator, levels=levels)


plt.colorbar(CS, label=r"R_trap", ticks=locator)
plt.show()
