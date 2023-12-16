import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy
import math 

# Aufgabe a
t, delta_t, Amp, delta_A = np.genfromtxt("test.txt", unpack=True)
t_Einheit = t * 50
t = unp.uarray(t_Einheit, delta_t)
Amp = unp.uarray(Amp, delta_A)

# print(f"t= {t}") # in mikrosek
# print(f"A= {Amp}") # in Volt

positiv_Amp = np.abs(unp.nominal_values(Amp))
fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(unp.nominal_values(t), unp.nominal_values(Amp), "x",label="Messwerte")
ax1.grid()
ax1.set_xlabel(r"$t\,[\mu$s]")
ax1.set_ylabel(r"$U\,$[V]")
ax1.legend(loc="best")
fig.savefig("plot_a.pdf")

# Aufgabe c
f, U_c, delta_U_C, T, delta_T = np.genfromtxt("DritteAufgabe.txt", unpack=True)
T_Einheit = T * 20  
U_c = unp.uarray(U_c, delta_U_C)
T = unp.uarray(T_Einheit, delta_T)
U_0 = 2.15

print(f"U_C= {U_c}") # in mikrosek
print(f"T= {T}") # in Volt

# x = np.linspace(0, 10, 1000)
# y = x ** np.sin(x)

# fig, (ax1, ax2) = plt.subplots(1, 2, layout="constrained")
# ax1.plot(x, y, label="Kurve")
# ax1.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
# ax1.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
# ax1.legend(loc="best")

# ax2.plot(x, y, label="Kurve")
# ax2.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
# ax2.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
# ax2.legend(loc="best")

# fig.savefig("build/plot.pdf")