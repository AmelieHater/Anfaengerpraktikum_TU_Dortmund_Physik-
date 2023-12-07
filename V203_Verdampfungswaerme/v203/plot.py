import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

#Bestimmung der Verdampfungswärme
T_WFinC , p_WFinmbar = np.genfromtxt("WoulffscheFlasche.txt", unpack=True)
p_aAinbar, T_aAinC = np.genfromtxt ("andereApparatur.txt", unpack=True)
p_0_in_bar = 985 * 10**-3 

#print("Temperatur Versuch 1 in celsius:",T_WFinC,"Druck Versuch 1 in millibar:", p_WFinmbar)
#print("Temperatur Versuch 2 in celsius:", T_aAinC, "Druck in Versuch 2 in bar:",p_aAinbar)
p_WFinbar_nom = p_WFinmbar * 1000
T_WFinK_nom = T_WFinC + 273.15
p_aAinbar_nom = p_aAinbar
T_aAinK_nom = T_aAinC + 273.15

T_WFinK_Fehler =  1 #Grad in Kelvin 
p_WFinbar_Fehler = 1 * 10**-3 #bar
p_aAinbar_Fehler = 0.5 #bar
T_aAinK_Fehler = 1 #Grad in Kelvin 

T_WFinK_Fehler_Array = np.ones((76,1))* T_WFinK_Fehler
T_aAinK_Fehler_Array = np.ones((14,1))* T_aAinK_Fehler
p_WFinbar_Fehler_Array = np.ones((76,1))* p_WFinbar_Fehler
p_aAinbar_Fehler_Array = np.ones((14,1))* p_aAinbar_Fehler

p_WFinbar = unp.uarray(p_WFinbar_nom, p_WFinbar_Fehler_Array)
p_aAinbar = unp.uarray(p_aAinbar_nom, p_aAinbar_Fehler_Array)
T_WFinK = unp.uarray(T_WFinK_nom, T_WFinK_Fehler_Array)
T_aAinK = unp.uarray(T_aAinK_nom, T_aAinK_Fehler_Array)

p_WFinbar_nom_durch_p_0 = p_WFinbar_nom / p_0_in_bar
#print("1.Versuch Druck in bar: ",p_WFinbar, "1. Versuch Temperatur in Kelvin", T_WFinK)
#print("2.Versuch Druck in bar: ", p_aAinbar,"2. Versuch Temperatur in Kelvin", T_aAinK)

#Plots 
#Plot für millibarMessung
x = (1/T_WFinK_nom)
y = np.log(p_WFinbar_nom_durch_p_0)

params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
#params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
#errors = np.sqrt(np.diag(covariance_matrix))

x_plot = np.linspace(0.0026,0.0034,10000)
plt.plot (x, y, "x", label = "Messwerte")
plt.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label ="Lineare Regression",
    linewidth=1,
)
plt.grid()
plt.xlabel(r'$T^{-1}$ [K $^{-1}$]')
plt.ylabel(r"$\ln{ \left ( \eta \right )}$")
plt.legend(loc = "best")
plt.savefig("plot1.pdf")