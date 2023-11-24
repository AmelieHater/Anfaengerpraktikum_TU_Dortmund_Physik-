import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

#Bestimmung Apperatekonstante D
Federwaagen_Abstand = 0.19975 #in m
phi_D, F_D = np.genfromtxt("D_Brchnng.txt", unpack=True)
#print(F_D * Federwaagen_Abstand)
D_0 = (F_D * Federwaagen_Abstand)/phi_D #in Nm
#print(D_0)
D_mean = np.mean(D_0)
#print(D_mean)
D_std = np.std(D_0)
D = ufloat(D_mean, D_std)
print("D ist: ",D)

#Eigentraegheitsmoment der Drillachse 
Abstand_a, T5_Periodendauer= np.genfromtxt("I_D_Brchnng.txt", unpack=True)
#print(Abstand_a)
a_quadrat = np.square(Abstand_a)
T_quadrat = np.square(T5_Periodendauer/5)
#print(a_quadrat)
#print(T_quadrat)


x = a_quadrat
y = T_quadrat

x_plot= np.linspace(0, 0.1)
params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
#errors = np.sqrt(np.diag(covariance_matrix))

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(x, y, "o", label="Messwerte")
ax.plot(x_plot, params[0] * x_plot + params[1],
       label="Lineare Regression")
ax.set_xlabel(r"Abstand$\,\,a^2 \,\,[m^2]$")
ax.set_ylabel(r"Periodendauer $\,\,T^2\,\,[s^2]$")
ax.legend(loc="best")
fig.savefig("plot.pdf")

#Lineare regression: y = ax+b a = params[0] b = params[1]
#Durchmesser der Gewichte berechnen mit Fehler
Blaues_Hoehe0, Blaues_Durchmesser0, Rotes_Hoehe0, Rotes_Durchmesser0 = np.genfromtxt("Gewichte_Maße.txt", unpack=True)
Gewichte_Durchmesser0 = (Blaues_Durchmesser0 + Rotes_Durchmesser0)/2
Gewichte_Durchmesser_mean = np.mean(Gewichte_Durchmesser0)

Gewichte_Durchmesser_std = np.std(Gewichte_Durchmesser0)
Gewichte_Durchmesser = ufloat(Gewichte_Durchmesser_mean, Gewichte_Durchmesser_std)
print("Durchmesser der Gewichte ist ",Gewichte_Durchmesser)
#Höhe der Gewichte ausrechnen mit Fehler 
Gewicht_Hoehe0 = (Blaues_Hoehe0 + Rotes_Hoehe0) /2
Gewicht_Hoehe_mean = np.mean(Gewicht_Hoehe0)
Gewicht_Hoehe_std = np.std(Gewicht_Hoehe0)
Gewicht_Hoehe = ufloat(Gewicht_Hoehe_mean, Gewicht_Hoehe_std)
print("Hoehe der Gewichte ist: ",Gewicht_Hoehe)
print("b ist: ", params[1])

I_D = (params[1]*D) / (4 * np.pi**2) - 2 * 0.26115 * np.square(Gewichte_Durchmesser/2) - (1/6) * 0.26115 * np.square(Gewicht_Hoehe)
Test1 = (params[1]*D) / (4 * np.pi**2)
Test2 = 2 * 0.26115 * np.square(Gewichte_Durchmesser/2) + (1/6) * 0.26115 * np.square(Gewicht_Hoehe)
print("I_D ist gleich: ", I_D)
print("Vorderer Summant ist: ", Test1, "Hinterer Summant ist: ", Test2)
#print(params[0],params[1])


#for name, value, error in zip ('ab', params, errors):
 #   print(f"{name} = {value:.3f} ± {error:.3f}")


#ax2.plot(x, y, label="Kurve")
#ax2.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
#ax2.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
#ax2.legend(loc="best")

