import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy
from scipy.optimize import curve_fit
import math 

saegezahn_f_kHz, saegezahn_f_kHz_Fehler, saegezahn_U_dB, saegezahn_U_dB_Fehler = np.genfromtxt("Messwerte_saegezahn.txt", unpack=True)
rechteck_f_kHz, rechteck_f_kHz_Fehler, rechteck_U_dB, rechteck_U_dB_Fehler = np.genfromtxt("Messwerte_rechteck.txt", unpack=True)
dreieck_f_kHz,dreieck_f_kHz_Fehler, dreieck_U_dB, dreieck_U_dB_Fehler = np.genfromtxt("Messwerte_dreieck.txt", unpack=True)

saegezahn_U_V_nom = 10 **(saegezahn_U_dB/20)
rechteck_U_V_nom = 10 **(rechteck_U_dB/20)
dreieck_U_V_nom = 10 **(dreieck_U_dB/20)

saegezahn_U_V_Fehler = 10 **(saegezahn_U_dB_Fehler/20)
rechteck_U_V_Fehler = 10 **(rechteck_U_dB_Fehler/20)
dreieck_U_V_Fehler = 10 **(dreieck_U_dB_Fehler/20)

saegezahn_U_V = unp.uarray(saegezahn_U_V_nom, saegezahn_U_V_Fehler)
rechteck_U_V = unp.uarray(rechteck_U_V_nom, rechteck_U_V_Fehler)
dreieck_U_V = unp.uarray(dreieck_U_V_nom, dreieck_U_V_Fehler)


def exponential_fct (x, a, b):
    return(a * x**(b))

popt_saeg, pcov_saeg = curve_fit(exponential_fct, saegezahn_f_kHz, saegezahn_U_V_nom)
a_saeg, b_saeg = popt_saeg
print(f"a_saeg: {a_saeg}, b_saeg: {b_saeg}")

x = np.linspace(0.01,120, 1000)
y_saeg = exponential_fct(x, a_saeg, b_saeg)

fig, ax = plt.subplots(1, 1)
ax.plot(saegezahn_f_kHz, saegezahn_U_V_nom, "x", label = "Messwerte")
ax.plot(x, y_saeg, label = "Fit" )
ax.grid()
ax.set_xlim([0,120])
ax.set_ylim([0,100])
ax.set_xlabel(r'$f$ [kHz]')
ax.set_ylabel(r"$U_S$ [V]")
ax.legend(loc = "best")
fig.savefig("plot1.pdf")

popt_recht, pcov_recht = curve_fit(exponential_fct, rechteck_f_kHz, rechteck_U_V_nom)
a_recht, b_recht = popt_recht
print(f"a_recht: {a_recht}, b_recht: {b_recht}")

x = np.linspace(0.01,400, 1000)
y_recht = exponential_fct(x, a_recht, b_recht)

fig2, ax2 = plt.subplots(1, 1)
ax2.plot(rechteck_f_kHz, rechteck_U_V_nom, "x", label = "Messwerte")
ax2.plot(x, y_recht, label = "Fit" )
ax2.grid()
ax2.set_xlim([0,400])
ax2.set_ylim([0,100])
ax2.set_xlabel(r'$f$ [kHz]')
ax2.set_ylabel(r"$U_R$ [V]")
ax2.legend(loc = "best")
fig2.savefig("plot2.pdf")

popt_drei, pcov_drei = curve_fit(exponential_fct, dreieck_f_kHz, dreieck_U_V_nom)
a_drei, b_drei = popt_drei
print(f"a_drei: {a_drei}, b_drei: {b_drei}")

x = np.linspace(0.01,140, 1000)
y_drei = exponential_fct(x, a_drei, b_drei)

fig3, ax3 = plt.subplots(1, 1)
ax3.plot(dreieck_f_kHz, dreieck_U_V_nom, "x", label = "Messwerte")
ax3.plot(x, y_drei, label = "Fit" )
ax3.grid()
ax3.set_xlim([0,140])
ax3.set_ylim([-10,70])
ax3.set_xlabel(r'$f$ [kHz]')
ax3.set_ylabel(r"$U_D$ [V]")
ax3.legend(loc = "best")
fig3.savefig("plot3.pdf")

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