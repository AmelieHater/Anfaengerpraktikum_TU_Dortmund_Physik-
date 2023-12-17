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
print(f"t pos.={t_pos} in mikrosek") 
print(f"A post= {Amp_pos} in V")

t_neg, delta_t, Amp_neg, delta_A = np.genfromtxt("negativeAmp.txt", unpack=True)
t_Einheit_neg =t_neg * 50
t_neg = unp.uarray(t_Einheit_neg+215, delta_t)
Amp_neg = unp.uarray(Amp_neg, delta_A)
print(f"t neg= {t_neg} in mikrosek")
print(f"A neg= {Amp_neg} in V")

# Plot für Aufgabe a 
def exponential_pos (x, a, b):
    return(a * np.exp(-x * b))

def exponential_neg (x, a, b):
    return(-a * np.exp(-x * b))

popt_pos, pcov_pos = curve_fit(exponential_pos, unp.nominal_values(t_pos), unp.nominal_values(Amp_pos))
a_pos, b_pos = popt_pos
print(f"a pos: {a_pos} in V")
print(f"b pos: {b_pos * 10**(6)} in Hz")

popt_neg, pcov_neg = curve_fit(exponential_neg, unp.nominal_values(t_neg), unp.nominal_values(Amp_neg))
a_neg, b_neg = popt_neg
print(f"a neg: {a_neg} in V")
print(f"b neg: {b_neg * 10**(6)} in Hz")

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

R_eff_exp = 2 * ((b_pos + b_neg)/2) * L * 10**6
print(f"R_eff_exp: {R_eff_exp} in ohm")

# ------------------- Aufgabe b -------------------
R_ap_exp = 4500
R_ap_theo = 2* unp.sqrt(L/C)
print(f"R_ap Theo: {R_ap_theo} (ohm)")
# ------------------- Aufgabe c -------------------
f, U_c, delta_U_C, T, delta_T = np.genfromtxt("DritteAufgabe.txt", unpack=True)
T_Einheit = T * 20  
U_c = unp.uarray(U_c, delta_U_C)
T = unp.uarray(T_Einheit, delta_T)
U_0 = 2.15
f_Einheit = f * 1000 # Einheit

U_zu_U0 = unp.nominal_values(U_c) /U_0
print(f"U/U_0: {U_zu_U0}")


def Kurve_c(omega, R, L, C):
    return (1/np.sqrt((1- L * C * omega**2)**2 + omega**2 * R**2 * C**2))

x = np.linspace(9.5, 45, 1000)
y_theoKurve = Kurve_c(2*np.pi*1000*x, unp.nominal_values(R2), unp.nominal_values(L), unp.nominal_values(C))
fig2, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(f, unp.nominal_values(U_zu_U0), "rx",label="Messwerte")
ax1.plot(x,y_theoKurve, "b-" , label = "Theoriekurve")
ax1.set_xscale("log")
ax1.grid(which="both")
ax1.set_xlim([9.5,45])
ax1.set_xlabel(r"$f\,\,$[kHz]")
ax1.set_ylabel(r"$\frac{U}{U_0}$")
ax1.legend(loc="best")
fig2.savefig("plot_b.pdf")
# print(f"U_C= {U_c} in V")
# print(f"T= {T} in mikrosek") 

print(f"Resonanzüberhöhung q theo: {max(y_theoKurve)}")

delta_omega = R2 / L 
delta_freq = delta_omega / (2 * np.pi)
omega_0 = unp.sqrt(1/(L*C))
q = omega_0/delta_omega
delta_freq_exp = 30.16 - 20.85 
print(f"delta_omega theo: {delta_omega} in Hz")
print(f'delta_freq  theo(Breite der Resonanzkurve): {delta_freq} in Hz')
print(f'delta_freq  exp(Breite der Resonanzkurve): {delta_freq_exp * 1000} in Hz')
print(f'Güte/Resonanzüberhöhung q theo: {q} in V')

U_omegaplus = max(y_theoKurve)/(np.sqrt(2))
print(f"Spannung von omega+ und omega- theo: {U_omegaplus}")
max_exp = ufloat(3.1395348837209305,0.023255813953488372)
U_omegaplus_exp = max_exp/(np.sqrt(2))
print(f"Spannung von omega+ und omega- exp: {U_omegaplus_exp}")

omega_res = unp.sqrt(1/(L*C) - (R2**2)/(2*(L**2)))
omega_1 = R2/(2*L) + unp.sqrt(R2**2/(4*(L**2)) + 1/(L*C))
omega_2 = -R2/(2*L) + unp.sqrt(R2**2/(4*(L**2)) + 1/(L*C))
print(f'Frequenz 1 theo: {omega_1/(1000*2*np.pi)} in kHz')
print(f'Frequenz 2 theo: {omega_2/(1000*2*np.pi)} in kHz')

def breite_Kurve(x, max):
    return(max/(np.sqrt(2))*x /x)
x1 = np.linspace(20.85, 30.16, 2)
fig3, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(f, unp.nominal_values(U_zu_U0), "r",label="Messwerte")
ax1.plot(x1, unp.nominal_values(breite_Kurve(x1,max_exp)), "b--" , label = "Experimentelle Breite")
ax1.grid(which="both")
ax1.set_xlim([18,33])
ax1.set_ylim([1,3.5])
ax1.set_xlabel(r"$f\,\,$[kHz]")
ax1.set_ylabel(r"$\frac{U}{U_0}$")
ax1.legend(loc="best")
fig3.savefig("plot_c.pdf")

# ------------------- Diskussion -------------------
def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo))
print("Diskussion:")
R_eff_abw = rel_Abweichung(R_eff_exp, R1)
print(f"rel. Abw. R eff: {R_eff_abw}")

R_ap_abw = rel_Abweichung(R_ap_exp, R_ap_theo)
print(f"rel. Abw. R ap: {R_ap_abw}")

max_abw = rel_Abweichung(unp.nominal_values(max_exp), max(y_theoKurve))
print(f"rel. Abw. Maximum(Resonanzüberhöhung q): {max_abw}")

breite_abw = rel_Abweichung(delta_freq_exp,delta_freq)
print(f"rel. Abw. Breite: {breite_abw}")

delta_freq_abw = rel_Abweichung(delta_freq_exp, delta_freq)
print(f"rel. Abw. delta_freq: {delta_freq_abw}")