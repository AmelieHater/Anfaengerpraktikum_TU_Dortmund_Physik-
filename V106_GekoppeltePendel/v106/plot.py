import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy

T_5_linkes, T_5_rechtes= np.genfromtxt("Laenge1_Einzelnes_Pendel_5T.txt", unpack=True)
T_linkes0 = T_5_linkes/5
T_linkes_mean = np.mean(T_linkes0)
T_linkes_std = scipy.stats.sem(T_linkes0)
T_linkes = ufloat(T_linkes_mean, T_linkes_std)
# print(f"T_linkes_L1: {T_linkes}")
T_rechtes0 = T_5_rechtes/5
T_rechtes_mean = np.mean(T_rechtes0)
T_rechtes_std = scipy.stats.sem(T_rechtes0)
T_rechtes = ufloat(T_rechtes_mean, T_rechtes_std)
# print(f"T_rechtes_L1: {T_rechtes}")

T_5_linkes_L2, T_5_rechtes_L2= np.genfromtxt("Laenge2_Einzelnes_Pendel_5T.txt", unpack=True)
T_linkes0_L2 = T_5_linkes_L2/5
T_linkes_L2_mean = np.mean(T_linkes0_L2)
T_linkes_L2_std = scipy.stats.sem(T_linkes0_L2)
T_linkes_L2 = ufloat(T_linkes_L2_mean, T_linkes_L2_std)
# print(f"T_linkes_L2: {T_linkes_L2}")
T_rechtes0_L2 = T_5_rechtes_L2/5
T_rechtes_L2_mean = np.mean(T_rechtes0_L2)
T_rechtes_L2_std = scipy.stats.sem(T_rechtes0_L2)
T_rechtes_L2 = ufloat(T_rechtes_L2_mean, T_rechtes_L2_std)
# print(f"T_rechtes_L2: {T_rechtes_L2}")

# Länge 1 gleichphasige Schwingung
g = 9.81
l_1 = 0.3000
omega_plus_1_theo = unp.sqrt(g/l_1)
print(f"omega_plus_1_theo: {omega_plus_1_theo}")
T_plus_1_theo = 2 * np.pi * unp.sqrt(l_1/g)
print(f"T_plus_1_theo: {T_plus_1_theo}")

# Länge 1 gegenphasige Schwingung


# Länge 2 gleichphasige Schwingung
l_2 = 0.6500
omega_plus_2_theo = unp.sqrt(g/l_2)
print(f"omega_plus_2_theo: {omega_plus_2_theo}")
T_plus_2_theo = 2 * np.pi * unp.sqrt(l_2/g)
print(f"T_plus_2_theo: {T_plus_2_theo}")

# Länge 2 gegenphasige Schwingung


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