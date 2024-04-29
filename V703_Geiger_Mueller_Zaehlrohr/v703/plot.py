import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
import scipy as scipy
import scipy.constants as const
from scipy.optimize import curve_fit

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent
def linear(x, m, b):
    return (m*x +b)

U_V, N_60s = np.genfromtxt("content/Messwerte/Kennlinie.txt", unpack = True)

N_s = N_60s / 60 # Zählrate pro sek
print (f"N_s: {N_s} 1/s")

U_A = 560 #V
x = np.linspace(380,700,10000) 

params, pcov = curve_fit(linear, U_V[4:21], N_s[4:21])
std = np.sqrt(np.diag(pcov))
m = unp.uarray(params[0], std[0])
b = unp.uarray(params[1], std[1])

print(f"m= {m}") # 1/Vs
print(f"b= {b}") # 1/s

fig1, ax1 = plt.subplots(1, 1, layout = "constrained")
ax1.plot(U_V, N_s, "x", label = "Messwerte")
ax1.plot(x, linear(x,params[0],params[1]), "-", label ="Lineare Regression")
ax1.plot(560, 16739/60, "o", fillstyle = "none", label = r"Arbeitspunkt $U_{\text{A}}$")
ax1.set_xlabel(r"$U\,$[V]")
ax1.set_ylabel(r"$N\,\left[\text{s}^{-1}\right]$")
# ax1.set_yscale("log")
ax1.legend(loc="best")
ax1.grid()
fig1.savefig("Kennlinie.pdf")

m_plateau = 100*100*m / ((16739/60))
z = 320 # Länge Plateau
s = ((U_A + 50) - (U_A - 50)) / U_A * 10
print(f"Plateusteigung s {s}")
print(f"Plateusteigung m {m_plateau}")

N_1 = (59355 + 67000 * 2) / 120 #1/s
N_2 = (67000 * 1 + 42837) / 120 
N_12 = (67000 * 3 + 65228) / 120 

print (f"N_1: {N_1 *120} 1/60s und {N_1} 1/s")
print (f"N_2: {N_2 *120} 1/60s und {N_2} 1/s")
print (f"N_12:{N_12 *120} 1/60s und {N_12} 1/s")

tau_2Q = (N_1 + N_2 - N_12) / (N_12**2 - N_1**2 - N_2**2) * 10**6
N_0 = N_12 /(1 - N_12 * tau_2Q)
verlust = N_0 - N_12
print(f"Totzeit: {tau_2Q } mu s")
print(f"tatsächliche Zählrate: {N_0} 1/s")
print(f"Zählratenverlust: {verlust} 1/s")

print(f"Diskussion:")
tau_rel = rel_Abweichung(225,tau_2Q)
print(f"rel. Abweichung tau: {tau_rel} %")

s_rel = rel_Abweichung(m_plateau,s)
print(f"rel. Abweichung Steigung s: {s_rel} %")
# fig2, ax2 = plt.subplots(1, 1, layout = "constrained")
# ax2.plot(x, linear(x,params[0],params[1]), "-", label ="Lineare Regression")
# ax2.plot(U_V, N_s, "x", label = "Messwerte")
# ax2.plot(560, 16739/60, "o", fillstyle = "none", label = r"Arbeitspunkt $U_{\text{A}}$")
# ax2.set_xlabel(r"$U\,$[V]")
# ax2.set_ylabel(r"$N\,\left[\text{s}^{-1}\right]$")
# ax2.legend(loc="best")
# ax2.grid()
# fig2.savefig("Linear.pdf")
