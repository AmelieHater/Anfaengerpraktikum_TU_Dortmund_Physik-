import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

#Bestimmung der Verdampfungswärme
T_WFinC , p_WFinmbar, T_WFinK_Fehler, p_WFinmbar_Fehler = np.genfromtxt("WoulffscheFlasche.txt", unpack=True)
p_aAinbar, T_aAinC, p_aAinbar_Fehler, T_aAinK_Fehler = np.genfromtxt ("andereApparatur.txt", unpack=True)
p_0_in_bar = 985 * 10**-3 

p_WFinbar_nom = p_WFinmbar *10**-3
T_WFinK_nom = T_WFinC + 273.15
p_aAinbar_nom = p_aAinbar
T_aAinK_nom = T_aAinC + 273.15

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
L_a = R * 373.15 # = 3102.541725963881 in Joul pro mol
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