import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.constants as const
from uncertainties import ufloat

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent


U_V_1, I_mA_1 = np.genfromtxt("content/Messdaten/Kennlinie_1.txt", unpack = True) 
U_V_2, I_mA_2 = np.genfromtxt("content/Messdaten/Kennlinie_2.txt", unpack = True) 
U_V_3, I_mA_3 = np.genfromtxt("content/Messdaten/Kennlinie_3.txt", unpack = True) 
U_V_4, I_mA_4 = np.genfromtxt("content/Messdaten/Kennlinie_4.txt", unpack = True) 
U_V_5, I_mA_5 = np.genfromtxt("content/Messdaten/Kennlinie_5.txt", unpack = True) 

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(U_V_1, I_mA_1, "x", label="bei $I_{Heiz} = 2$")
ax1.plot(U_V_2, I_mA_2, "x", label="bei $I_{Heiz} = 2,1$")
ax1.plot(U_V_3, I_mA_3, "x", label="bei $I_{Heiz} = 2,2$")
ax1.plot(U_V_4, I_mA_4, "x", label="bei $I_{Heiz} = 2,3$")
ax1.set_xlabel(r"$U \,\,[\text{V}]$")
ax1.set_ylabel(r"$A \,\,[\text{mA}]$")
ax1.legend(loc="best")

fig.savefig("plot_1.pdf")


fig, (ax2) = plt.subplots(1, 1, layout="constrained")
ax2.plot(U_V_5, I_mA_5, "x",  label="bei $I_{Heiz} = 2,4$")
ax2.set_xlabel(r"$U \,\,[\text{V}]$")
ax2.set_ylabel(r"$A \,\,[\text{mA}]$")
ax2.legend(loc="best")

fig.savefig("plot_2.pdf")

#Sättigungsstrom: 
#1.Kennlinie = 0.060
#2.Kennlinie = 0.140
#3.Kennlinie = 0.319
#4.Kennlinie = 0.719
#5.Kennlinie = 1.356

#Gültigkeit der Langmuir-Schottkyschen Gleichung
#Regression
# Fit a polynomial of degree 1, return covariance matrix
params_1, covariance_matrix_1 = np.polyfit(np.log(U_V_5), np.log(I_mA_5), deg=1, cov=True)

errors_1 = np.sqrt(np.diag(covariance_matrix_1))

for name, value, error in zip("ab", params_1, errors_1):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x = np.linspace(1,5.5,10)
fig, (ax3) = plt.subplots(1, 1, layout="constrained")
ax3.plot(np.log(U_V_5), np.log(I_mA_5), "x",  label="Messwerte")
ax3.plot(x, params_1[0]*x + params_1[1], label="Lineare Regression" )
ax3.set_xlabel(r"$\text{log}(U / \text{V})$")
ax3.set_ylabel(r"$\text{log}(I / \text{mA})$")
ax3.legend(loc="best")

fig.savefig("plot_3.pdf")

#Als a und b ergeben sich diese Werte:
#y = ax + b
#a = 1.256 ± 0.020
#b = -6.277 ± 0.078

#Bestimmung der Kathodentemperatur durch Untersuchung des Anlaufstromgebiets
U_gegen_V, I_nA = np.genfromtxt("content/Messdaten/Anlaufkurve.txt", unpack = True)
params_2, covariance_matrix_2 = np.polyfit(U_gegen_V, np.log(I_nA), deg=1, cov=True)

errors_2 = np.sqrt(np.diag(covariance_matrix_2))

for name, value, error in zip("cd", params_2, errors_2):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_2 = np.linspace(-0.01,0.66,10)
fig, (ax4) = plt.subplots(1, 1, layout="constrained")
ax4.plot(U_gegen_V, np.log(I_nA), "x",  label="Messwerte")
ax4.plot(x_2, params_2[0]*x_2 + params_2[1], label="Lineare Regression" )
ax4.set_xlabel(r"$U \,\,[\text{V}]$")
ax4.set_ylabel(r"$\text{log}(I / \text{nA})$")
ax4.legend(loc="best")

fig.savefig("plot_4.pdf")

#Lineare Regression ergibt: mit y = c *x + d
#c = -7.819 ± 0.626
#d = 1.765 ± 0.225

#Kathodentemeratur ergibt sich dann zu t = - e/(kc)
c = ufloat(params_2[0],errors_2[0])
T = - const.elementary_charge / (const.Boltzmann * c)
print("T = ",T)
#T =  1484.179499962159

#Bestimmung der Kathodentemperaturen und der Austrittsarbeit von Wolfram
sigma = 5.7 * 10**(-12) # in W/(cm² K⁴)
f = 0.32 #in cm²
eta = 0.28 #Emissionsgrad der Oberfläche 
N_WL = 0.95 #in Watt
I_Heiz = (2,2.1,2.2,2.3,2.4)
U_Heiz = (4,4,4.5,5,5)
T_array = ((np.multiply(I_Heiz, U_Heiz) - N_WL)/(f*eta*sigma))**(1/4)
print(T_array)
#Output: [1927.53184987 1954.30953913 2046.0207907  2131.90279307 2156.72540037]

#Austrittsarbeit damit berechnen
#1.Kennlinie = 0.060
#2.Kennlinie = 0.140
#3.Kennlinie = 0.319
#4.Kennlinie = 0.719
#5.Kennlinie = 1.356
I_S = np.array([0.060, 0.140, 0.319, 0.719, 1.356]) #ist wahrscheinlich der Sättigungsstrom
phi = - (const.Boltzmann * T_array)/(const.elementary_charge) * np.log(I_S * (const.Planck)**3 /(4 * np.pi * f * const.electron_mass * const.elementary_charge * (const.Boltzmann)**2 * T_array**2))
print(const.Planck)
print(phi)
#Output phi = [5.11614313 5.04917175 5.15708723 5.23936631 5.18676248]
#Mittelwert der Austrittsarbeiten bestimmen 
phi_mean = ufloat(np.mean(phi), sp.stats.sem(phi))
print(phi_mean) #Output: 5.150+/-0.032

#Relative Abweichung bestimmen
print("Abweichung zu 3/2: ", rel_Abweichung(1.26, 3/2))
print("Abweichung der Kathodentemperatur: ", rel_Abweichung(1500, 2131.90))
print("Abweichung der Austrittsarbeit: ", rel_Abweichung(5.15, 4.54))

#Abweichung zu 3/2:  16.0
#Abweichung der Kathodentemperatur:  29.640227027534127
#Abweichung der Austrittsarbeit:  13.43612334801763