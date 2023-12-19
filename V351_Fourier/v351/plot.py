import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import math 

saegezahn_f_kHz, saegezahn_U_dB = np.genfromtxt("Messwerte_saegezahn.txt", unpack=True)
rechteck_f_kHz, rechteck_U_dB = np.genfromtxt("Messwerte_rechteck.txt", unpack=True)
dreieck_f_kHz, dreieck_U_dB = np.genfromtxt("Messwerte_dreieck.txt", unpack=True)

saegezahn_U_V = 10 **(saegezahn_U_dB/20)
rechteck_U_V = 10 **(rechteck_U_dB/20)
dreieck_U_V = 10 **(dreieck_U_dB/20)


x = np.linspace(0, 10, 1000)
y = x ** np.sin(x)
fig, ax = plt.subplots(1, 1)
ax.plot(x, y, "x", label = "Messwerte")
ax.grid()
ax.set_xlim([10,110])
ax.set_xlabel(r'$f$ [kHz]')
ax.set_ylabel(r"$U$ [V]")
ax.legend(loc = "best")
fig.savefig("plot.pdf")

print("Saegezahn", saegezahn_f_kHz)
print(saegezahn_U_V)
#x = np.linspace(0, 10, 1000)
#y = x ** np.sin(x)
#
#fig, (ax1, ax2) = plt.subplots(1, 2, layout="constrained")
#ax1.plot(x, y, label="Kurve")
#ax1.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
#ax1.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
#ax1.legend(loc="best")
#
#ax2.plot(x, y, label="Kurve")
#ax2.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
#ax2.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
#ax2.legend(loc="best")
#
#fig.savefig("build/plot.pdf")