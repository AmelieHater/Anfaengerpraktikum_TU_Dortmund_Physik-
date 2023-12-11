import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import math 

#Bestimmung der Verdampfungswärme
T_WFinC , p_WFinmbar, T_WFinK_Fehler, p_WFinmbar_Fehler = np.genfromtxt("WoulffscheFlasche.txt", unpack=True)
p_0_in_bar = 985 * 10**-3 

p_WFinbar_nom = p_WFinmbar *10**-3
T_WFinK_nom = T_WFinC + 273.15


p_WFinbar_Fehler = p_WFinmbar_Fehler *10**-3


#T_WFinK_Fehler =  1 #Grad in Kelvin 
#p_WFinbar_Fehler = 1 * 10**-3 #bar
#p_aAinbar_Fehler = 0.5 #bar
#T_aAinK_Fehler = 1 #Grad in Kelvin 

T_WFinK = unp.uarray(T_WFinK_nom, T_WFinK_Fehler)
p_WFinbar = unp.uarray(p_WFinbar_nom, p_WFinbar_Fehler)
#print("Temperatur: ",T_WFinK)
#print("Druck: ", p_WFinbar)


p_WFinbar_durch_p_0 = p_WFinbar / p_0_in_bar


#Plot für millibarMessung
x = (1/T_WFinK)
y = unp.log(p_WFinbar_durch_p_0)

#print("y Werte: ",y)
#print(T_WFinK_nom)
#print(1/T_WFinK_nom)

params, covariance_matrix = np.polyfit(unp.nominal_values(x), unp.nominal_values(y), deg=1, cov=True)
errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip ('mn', params, errors):
    print(f"{name} = {value:.3f} ± {error:.3f}")

fig, ax = plt.subplots(1, 1)
x_plot = np.linspace(0.00265,0.00337,10000)
ax.plot (unp.nominal_values(x), unp.nominal_values(y), "x", label = "Messwerte")
ax.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label ="Lineare Regression",
    linewidth=1,
)
ax.grid()
ax.set_xlim([0.00265,0.00337])
ax.set_xlabel(r'$T^{-1}$ [K $^{-1}$]')
ax.set_ylabel(r"$\ln{ \left ( \frac{p}{p_0} \right )}$")
ax.legend(loc = "best")
fig.savefig("plot1.pdf")

#Ausgleichrechnung zur Bestimmung von L y=mx+n
m = ufloat(-3416.886 , 52.501) # m = -3416.8862761957876
n = ufloat(8.997, 0.157)

#Berechnung L
R = 8.31446261815324 # in Joul pro(mol * Kelvin)
L = -m * R # = (2.84+/-0.04)e+04 in Joul pro mol
#print(L) 

#Berechnung L_a = R * T (L_a äußere Verdampfungswärme bei T = 100 °C)
L_a = R * 373 # = 3101.2945565711584in Joul pro mol
#print(L_a) 

#Berechnung von L_i = L - L_a 
L_i = L - L_a # (2.53+/-0.04)e+04 in Joul pro mol
#print(L_i)

#Berechnung von L_i pro Molekül, dh. durch Avogadrokonstante N_A teilen 
N_A = 6.02214076 * 10 **23 # in 1/mol
L_i_M_Joule = L_i/N_A # = (4.20+/-0.07)e-20 in Joul
#print(L_i_M_Joule)
Elementarladung = 1.602176634 * 10 **-19 # in As
L_i_M_eV = L_i_M_Joule /Elementarladung # 0.262+/-0.005 in eV
#print(L_i_M_eV)


# 2. Versuch: Temperaturabhängigkeit der Verdampfungswärme 
p_aAinbar, T_aAinC, p_aAinbar_Fehler, T_aAinK_Fehler = np.genfromtxt ("andereApparatur.txt", unpack=True)

p_aAinkilopascal_nom = p_aAinbar * 100
T_aAinK_nom = T_aAinC + 273.15
p_aAinkilopascal_Fehler = p_aAinbar_Fehler * 100

T_aAinK = unp.uarray(T_aAinK_nom, T_aAinK_Fehler)
p_aAinkilopascal = unp.uarray(p_aAinkilopascal_nom, p_aAinkilopascal_Fehler)

#Plot für millibarMessung
x = T_aAinK
y = p_aAinkilopascal


params, covariance_matrix = np.polyfit(unp.nominal_values(x), unp.nominal_values(y), deg=3, cov=True)
errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip ('abcd', params, errors):
    print(f"{name} = {value:.6f} ± {error:.6f}")

fig2, ax = plt.subplots(1, 1)
x_plot = np.linspace(375,475,10000)
ax.plot (unp.nominal_values(x), unp.nominal_values(y), "x", label = "Messwerte")
ax.plot(
    x_plot,
    params[0] * (x_plot)**3 + params[1] * (x_plot)**2 + params[2] * x_plot + params[3],
    label = "Regression",
    linewidth=1,
)
ax.grid()
ax.set_xlim([375,475])
ax.set_xlabel(r'$T$ [K]')
ax.set_ylabel(r"$p$ [kPa]")
ax.legend(loc = "best")
fig2.savefig("plot2.pdf")

# Polynom: f(x) = ax³ + bx² + cx + d
a_mit_Fehler = ufloat(0.000644, 0.000211)
b_mit_Fehler = ufloat(-0.661278, 0.273298)
c_mit_Fehler = ufloat(227.259041, 117.624126)
d_mit_Fehler = ufloat(-26160.860147, 16842.766350)

a = 0.000644
b = -0.661278
c = 227.259041
d = -26160.860147

#Nochmal 2 riesige Plots yey 
#Zuerst für L_-(T) von V_-
a_Schlange = 0.9 # Jm³ / mol²

fig3, ax = plt.subplots(1, 1)
x_plot = np.linspace(375,475,10000)
p = a * x_plot**3 + b * x_plot**2 + c * x_plot + d
p_Ableitung = 3*a * x_plot**2 + 2 * b * x_plot + c

y_plot_1 = x_plot * (((R * x_plot)/(2*p)) - np.sqrt(((R * x_plot)/(2*p))**2 - (a_Schlange/p))) * p_Ableitung

ax.plot(
    x_plot,
    y_plot_1,
    linewidth=1,
)
ax.grid()
ax.set_xlim([375,475])
ax.set_xlabel(r'$T$ [K]')
ax.set_ylabel(r"$L_-(T)$ [Jmol$^{-1}]$")
fig3.savefig("plot3.pdf")

#Dann für L_+(T) von V_+

fig4, ax = plt.subplots(1, 1)
x_plot = np.linspace(375,475,10000)


y_plot_2 = x_plot * (((R * x_plot)/(2*p)) + np.sqrt(((R * x_plot)/(2*p))**2 - (a_Schlange/p))) * p_Ableitung 
ax.plot(
    x_plot,
    y_plot_2 * 10 **-3,
    linewidth=1,
)
ax.grid()
ax.set_xlim([375,475])
ax.set_xlabel(r'$T$ [K]')
ax.set_ylabel(r"$L_+(T)$ [$\cdot$ 10$^{3}$ Jmol$^{-1}$]")
fig4.savefig("plot4.pdf")