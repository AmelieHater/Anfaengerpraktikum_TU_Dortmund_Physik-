import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy
from scipy.optimize import curve_fit
import math 

# ------------------- Aufgabe a -------------------
L = ufloat(16.87, 0.05) * 10**(-3)
R1 = ufloat(67.2, 0.1)
R2 = ufloat(682, 0.5)
C = ufloat(2.060, 0.003) * 10**(-9)
t_pos, delta_t, Amp_pos, delta_A = np.genfromtxt("positiveAmp.txt", unpack=True)
t_Einheit_pos =t_pos * 50
t_pos = unp.uarray(t_Einheit_pos+215, delta_t)
Amp_pos = unp.uarray(Amp_pos, delta_A)
# print(f"t pos.={t_pos}") # in mikrosek
# print(f"A post= {Amp_pos}") # in Volt 

t_neg, delta_t, Amp_neg, delta_A = np.genfromtxt("negativeAmp.txt", unpack=True)
t_Einheit_neg =t_neg * 50
t_neg = unp.uarray(t_Einheit_neg+215, delta_t)
Amp_neg = unp.uarray(Amp_neg, delta_A)
# print(f"t neg= {t_neg}") # in mikrosek
# print(f"A neg= {Amp_neg}") # in Volt

# Plot für Aufgabe a 
def exponential_pos (x, a, b):
    return(a * np.exp(-x * b))

def exponential_neg (x, a, b):
    return(-a * np.exp(-x * b))

popt_pos, pcov_pos = curve_fit(exponential_pos, unp.nominal_values(t_pos), unp.nominal_values(Amp_pos))
a_pos, b_pos = popt_pos
print(f"a: {a_pos}")
print(f"b: {b_pos}")

popt_neg, pcov_neg = curve_fit(exponential_neg, unp.nominal_values(t_neg), unp.nominal_values(Amp_neg))
a_neg, b_neg = popt_neg
print(f"a: {a_neg}")
print(f"b: {b_neg}")

x = np.linspace(0,430, 100000)
y_pos = exponential_pos(x, a_pos, b_pos)
y_neg =exponential_neg(x,a_neg, b_neg)
fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(unp.nominal_values(t_pos), unp.nominal_values(Amp_pos), "rx",label="Messwerte")
ax1.plot(unp.nominal_values(t_neg), unp.nominal_values(Amp_neg), "rx")
ax1.plot(x, y_pos, "b-", label="Exponentieller Fit")
ax1.plot(x, y_neg, "g-", label="Exponentieller Fit")
ax1.grid()
ax1.set_xlim([0,430])
ax1.set_xlabel(r"$t\,\,[10^{-6}\,$s]")
ax1.set_ylabel(r"$U\,\,$[V]")
ax1.legend(loc="best")
fig.savefig("plot_a.pdf")

R_eff_exp = 2 * ((b_pos + b_neg)/2) * L * 10**3
print(f"R_eff_exp: {R_eff_exp}")

# ------------------- Aufgabe b -------------------
R_ap_theo = 2* unp.sqrt(L/C)
print(f"R_ap Theo: {R_ap_theo} (ohm)")
# ------------------- Aufgabe c -------------------
f, U_c, delta_U_C, T, delta_T = np.genfromtxt("DritteAufgabe.txt", unpack=True)
T_Einheit = T * 20  
U_c = unp.uarray(U_c, delta_U_C)
T = unp.uarray(T_Einheit, delta_T)
U_0 = 2.15
f_Einheit = f * 1000 # Einheit

U_zu_U0 = U_c /U_0
print(f"U/U_0: {U_zu_U0}")


def Kurve_c(omega, R, L, C):
    return (1/np.sqrt((1- L * C * omega**2)**2 + omega**2 * R**2 * C**2))

x = np.linspace(9.5, 45, 1000)
fig2, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(f, unp.nominal_values(U_zu_U0), "rx",label="Messwerte")
ax1.plot(x, Kurve_c(2*np.pi*1000*x, unp.nominal_values(R2), unp.nominal_values(L), unp.nominal_values(C)), "b-" , label = "Theoriekurve")
ax1.set_xscale("log")
ax1.grid(which="both")
ax1.set_xlim([9.5,45])
ax1.set_xlabel(r"$f\,\,$[kHz]]")
ax1.set_ylabel(r"$U/U_0$")
ax1.legend(loc="best")
fig2.savefig("plot_b.pdf")
# print(f"U_C= {U_c}") # in mikrosek
# print(f"T= {T}") # in Volt

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