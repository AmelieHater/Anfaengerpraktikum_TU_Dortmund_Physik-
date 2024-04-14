import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy
import scipy.constants as const
from scipy.optimize import curve_fit
# Konstanen und Funktionen
R = 13.9 #in eV (Rydbergenergie)
d = 201.4 * 10**(-12) # in m Gitterkonsttante
h = const.h
e = const.e
c = const.c
alpha = 7.287*10**(-3)
# Wellenlänge
def lambda_min(theta):
    return 2*d * np.sin(theta*np.pi/180)
# Abschimrkonstante
def sigma_k (Z,E):
    return Z - np.sqrt((E/R) - (((alpha**2)*(Z**4))/4))
# Energie
def energie_k(lambda_m):
    return h * c /(lambda_m * e)

# Bragg-Bedingung
theta2_Bragg, r_Bragg = np.genfromtxt("content/Messwerte/Bragg.txt", unpack =True)
theta_Bragg = theta2_Bragg / 2

plt.plot(theta_Bragg, r_Bragg, "x-", label="Messwerte")
plt.plot(14, 191, "o", color="red",label="Maximum")
plt.axvline(x=14,color='green' ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([13,15])
plt.ylim([40,200])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Bragg.pdf")
plt.close()

# Emissionsspektrum
theta2_Emission, r_Emission = np.genfromtxt("content/Messwerte/Emissionsspektrum.txt", unpack=True)
theta_Emission = theta2_Emission / 2

plt.plot(theta_Emission, r_Emission, "x", label="Messwerte")
plt.plot(theta_Emission[8:79], r_Emission[8:79], "-", color="orange", label=r"Bremsberg")
plt.plot(theta_Emission[90:96], r_Emission[90:96], "-", color="red", label=r"$K_{\alpha}$-Linie")
plt.plot(theta_Emission[79:84], r_Emission[79:84], "-", color="green", label=r"$K_{\beta}$-Linie")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([4,26])
plt.ylim([0,4900])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Emissionsspektrum.pdf")
plt.close()

lambda_min1 = lambda_min(theta_Emission[8])
energie_max = energie_k(lambda_min1)
print ("Auswertung Emissionssektrum:")
print (f"Grenzwinkel: {theta_Emission[8]}°")
print (f"min. Wellenänge: {lambda_min1} m")
print (f"max. Energie Bremsberg: {energie_max} eV")

# Detailspektrum
plt.plot(theta_Emission[70:100], r_Emission[70:100], "x-", label="Messwerte")
# plt.plot(theta_Emission[90:96], r_Emission[90:96], "-", color="red", label=r"$K_{\alpha}$-Linie")
# plt.plot(theta_Emission[79:84], r_Emission[79:84], "-", color="green", label=r"$K_{\beta}$-Linie")
plt.hlines(y=4756/2, xmin=22.30, xmax=22.75,linestyle = ":", color="red", label=r"$\text{Halbwertsbreite des } K_{\alpha}\text{-Peaks}$")
plt.hlines(y=1387/2, xmin = 20.10, xmax = 20.50, linestyle = ":", color="green", label=r"$\text{Halbwertsbreite des } K_{\beta}\text{-Peaks}$")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([19,24])
plt.ylim([0,4900])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Detailspektrum.pdf")
plt.close()

lambda_K_beta = lambda_min(20.2)
lambda_K_beta1 = lambda_min(20.1)
lambda_K_beta2 = lambda_min(20.5)
Energie_K_beta1 = energie_k(lambda_K_beta1)
Energie_K_beta2 = energie_k(lambda_K_beta2)
Energie_K_beta = energie_k(lambda_K_beta)

lambda_K_alpha = lambda_min(22.6)
lambda_K_alpha1 = lambda_min(22.3)
lambda_K_alpha2 = lambda_min(22.75)
Energie_K_alpha1 = energie_k(lambda_K_alpha1)
Energie_K_alpha2 = energie_k(lambda_K_alpha2)
Energie_K_alpha = energie_k(lambda_K_alpha)

print("Auswertung Detailspektrum")
print(f"K_beta Linie Halbwertsbreite: 0,4°")
print(f"K_beta Linie min: 20,10° mit Energie {Energie_K_beta1}eV")
print(f"K_beta Linie max: 20,50° mit Energie {Energie_K_beta2}eV")
print(f"K_beta Linie max: 20,20° mit Energie {Energie_K_beta}eV")
print(f"K_beta Auflösungsvermögen A: {Energie_K_beta/(Energie_K_beta1 - Energie_K_beta2)}")
print ("----------------")
print(f"K_alpha Linie Halbwertsbreite: 0,45°")
print(f"K_alpha Linie min: 22,30° mit Energie {Energie_K_alpha1}eV")
print(f"K_alpha Linie max: 22,75° mit Energie {Energie_K_alpha2}eV")
print(f"K_alpha Linie max: 22,60° mit Energie {Energie_K_alpha}eV")
print(f"K_alpha Auflösungsvermögen A: {Energie_K_alpha/(Energie_K_alpha1 - Energie_K_alpha2)}")

# Abschirmkonstanten
sigma_1_exp = 29 - np.sqrt(8.988*1000/R)
sigma_2_exp = 29 - np.sqrt(4*((29 - sigma_1_exp)**2 - (Energie_K_alpha/R)))
sigma_3_exp = 29 - np.sqrt(9*((29 - sigma_1_exp)**2 - (Energie_K_beta/R)))
print("----------------")
print(f"sigma_1 exp: {sigma_1_exp}")
print(f"sigma_2 exp: {sigma_2_exp}")
print(f"sigma_3 exp: {sigma_3_exp}")
print("----------------")
sigma_1_theo = sigma_1_exp
sigma_2_theo = 29 - np.sqrt(4*((29 - sigma_1_exp)**2 - (8 *1000/R)))
sigma_3_theo = 29 - np.sqrt(9*((29 - sigma_1_exp)**2 - (8.95*1000/R)))
print(f"sigma_1 theo: {sigma_1_theo}")
print(f"sigma_2 theo: {sigma_2_theo}")
print(f"sigma_3 theo: {sigma_3_theo}")

# Strontium
theta2_Sr, r_Sr = np.genfromtxt("content/Messwerte/Strontium.txt", unpack=True)
theta_Sr = theta2_Sr / 2

mittel_unten_Sr = np.mean(r_Sr[0:8])
mittel_oben_Sr = np.mean(r_Sr[12:20])
mittel_Sr = np.mean(r_Sr)

plt.plot(theta_Sr, r_Sr, "x-", label="Messwerte")
plt.hlines(y=mittel_unten_Sr, xmin = theta_Sr[0], xmax = theta_Sr[8], linestyle=":", color = "red" , label = r"Mittelwert unten")
plt.hlines(y=mittel_oben_Sr, xmin = theta_Sr[12], xmax = theta_Sr[20], linestyle=":", color = "orange" , label = r"Mittelwert oben")
plt.hlines(y=mittel_Sr, xmin = theta_Sr[0], xmax = theta_Sr[20], linestyle="--", color = "green" , label = r"Mittelwert")
plt.axvline(x = 10.99, linestyle = "--", color = "purple", label = r"$\theta_{K} = 10,99\,°$" )
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([10,12])
plt.ylim([30,120])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Strontium.pdf")
plt.close()

Theta_K_Sr_exp = 10.99
print("----------------")
print(f"Theta_K_Sr_exp: {Theta_K_Sr_exp}°")

# Zirkonium
theta2_Zr, r_Zr = np.genfromtxt("content/Messwerte/Zirkonium.txt", unpack=True)
theta_Zr = theta2_Zr / 2

mittel_unten_Zr = np.mean(r_Zr[0:6])
mittel_oben_Zr = np.mean(r_Zr[12:19])
mittel_Zr = np.mean(r_Zr)

plt.plot(theta_Zr, r_Zr, "x-", label="Messwerte")
plt.hlines(y=mittel_unten_Zr, xmin = theta_Zr[0], xmax = theta_Zr[6], linestyle=":", color = "red" , label = r"Mittelwert unten")
plt.hlines(y=mittel_oben_Zr, xmin = theta_Zr[12], xmax = theta_Zr[19], linestyle=":", color = "orange" , label = r"Mittelwert oben")
plt.hlines(y=mittel_Zr, xmin = theta_Zr[0], xmax = theta_Zr[19], linestyle="--", color = "green" , label = r"Mittelwert")
plt.axvline(x = 9.91, linestyle = "--", color = "purple", label = r"$\theta_{K} = 9,91\,°$" )
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([9,10.9])
plt.ylim([70,180])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Zirkonium.pdf")
plt.close()

Theta_K_Zr_exp = 9.91
print(f"Theta_K_Zr_exp: {Theta_K_Zr_exp}°")

# Brom
theta2_Br, r_Br = np.genfromtxt("content/Messwerte/Brom.txt", unpack=True)
theta_Br = theta2_Br / 2

mittel_unten_Br = np.mean(r_Br[0:10])
mittel_oben_Br = np.mean(r_Br[13:19])
mittel_Br = np.mean(r_Br)

plt.plot(theta_Br, r_Br, "x-", label="Messwerte")
plt.hlines(y=mittel_unten_Br, xmin = theta_Br[0], xmax = theta_Br[10], linestyle=":", color = "red" , label = r"Mittelwert unten")
plt.hlines(y=mittel_oben_Br, xmin = theta_Br[13], xmax = theta_Br[19], linestyle=":", color = "orange" , label = r"Mittelwert oben")
plt.hlines(y=mittel_Br, xmin = theta_Br[0], xmax = theta_Br[19], linestyle="--", color = "green" , label = r"Mittelwert")
plt.axvline(x = 13.13, linestyle = "--", color = "purple", label = r"$\theta_{K} = 13,13\,°$" )
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([12,13.9])
plt.ylim([12,30])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Brom.pdf")
plt.close()

Theta_K_Br_exp = 13.13
print(f"Theta_K_Br_exp: {Theta_K_Br_exp}°")

# Zink
theta2_Zn, r_Zn = np.genfromtxt("content/Messwerte/Zink.txt", unpack=True)
theta_Zn = theta2_Zn / 2

# mittel_unten_Zn = np.mean(r_Zn[0:10])
# mittel_oben_Zn = np.mean(r_Zn[13:15])
mittel_Zn = np.mean(r_Zn)


plt.plot(theta_Zn, r_Zn, "x-", label="Messwerte")
# plt.hlines(y=mittel_unten_Zn, xmin = theta_Zn[0], xmax = theta_Zn[10], linestyle=":", color = "red" , label = r"Mittelwert unten")
# plt.hlines(y=mittel_oben_Zn, xmin = theta_Zn[13], xmax = theta_Zn[15], linestyle=":", color = "orange" , label = r"Mittelwert oben")
plt.hlines(y=mittel_Zn, xmin = theta_Zn[0], xmax = theta_Zn[15], linestyle="--", color = "green" , label = r"Mittelwert")
plt.axvline(x = 19.84, linestyle = "--", color = "purple", label = r"$\theta_{K} = 19,84\,°$" )
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([19.4,20.9])
plt.ylim([25,430])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Zink.pdf")
plt.close()

Theta_K_Zn_exp = 19.84
print(f"Theta_K_Zn_exp: {Theta_K_Zn_exp}°")
# Gallium
theta2_Ga, r_Ga = np.genfromtxt("content/Messwerte/Gallium.txt", unpack=True)
theta_Ga = theta2_Ga / 2

# mittel_unten_Ga = np.mean(r_Ga[0:10])
# mittel_oben_Ga = np.mean(r_Ga[13:20])
mittel_Ga = np.mean(r_Ga)

plt.plot(theta_Ga, r_Ga, "x-", label="Messwerte")
# plt.hlines(y=mittel_unten_Ga, xmin = theta_Ga[0], xmax = theta_Ga[10], linestyle=":", color = "red" , label = r"Mittelwert unten")
# plt.hlines(y=mittel_oben_Ga, xmin = theta_Ga[13], xmax = theta_Ga[20], linestyle=":", color = "orange" , label = r"Mittelwert oben")
plt.hlines(y=mittel_Ga, xmin = theta_Ga[0], xmax = theta_Ga[20], linestyle="--", color = "green" , label = r"Mittelwert")
plt.axvline(x = 17.21, linestyle = "--", color = "purple", label = r"$\theta_{K} = 17,21\,°$" )
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([17,19])
plt.ylim([20,34])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Gallium.pdf")
plt.close()

Theta_K_Ga_exp = 17.21
print(f"Theta_K_Ga_exp: {Theta_K_Ga_exp}°")
print("----------------")

# Energien  
energie_Sr_exp = energie_k(lambda_min(Theta_K_Sr_exp)) #eV
energie_Zr_exp = energie_k(lambda_min(Theta_K_Zr_exp)) #eV
energie_Br_exp = energie_k(lambda_min(Theta_K_Br_exp))#eV
energie_Zn_exp = energie_k(lambda_min(Theta_K_Zn_exp))#eV
energie_Ga_exp = energie_k(lambda_min(Theta_K_Ga_exp))#eV

print(f"E_Sr exp: {energie_Sr_exp} eV")
print(f"E_Zr exp: {energie_Zr_exp} eV")
print(f"E_Br exp: {energie_Br_exp} eV")
print(f"E_Zn exp: {energie_Zn_exp} eV")
print(f"E_Ga exp: {energie_Ga_exp} eV")

# Abschirmkonstante 
sigma_Sr_exp = sigma_k(38,energie_Sr_exp)
sigma_Zr_exp = sigma_k(40,energie_Zr_exp)
sigma_Br_exp = sigma_k(35,energie_Br_exp)
sigma_Zn_exp = sigma_k(30,energie_Zn_exp)
sigma_Ga_exp = sigma_k(31,energie_Ga_exp)

print("----------------")
print(f"sigma_Sr exp: {sigma_Sr_exp}")
print(f"sigma_Zr exp: {sigma_Zr_exp}")
print(f"sigma_Br exp: {sigma_Br_exp}")
print(f"sigma_Zn exp: {sigma_Zn_exp}")
print(f"sigma_Ga exp: {sigma_Ga_exp}")

# Rydbergkonstante
energie = [energie_Zn_exp, energie_Ga_exp, energie_Br_exp, energie_Sr_exp, energie_Zr_exp]
energie_sqrt = np.sqrt(energie)

Z = [30-sigma_Zn_exp, 31-sigma_Ga_exp, 35-sigma_Br_exp, 38-sigma_Sr_exp, 40-sigma_Zr_exp]

def linear(x,m,b):
    return m*x+b

params, pcov = curve_fit(linear, Z, energie_sqrt)
std = np.sqrt(np.diag(pcov))

a_sqrt = ufloat(params[0], std[0])
b_sqrt = ufloat(params[1], std[1])
a = a_sqrt**2
b = b_sqrt**2

print("----------------")
print(f"a_sqrt: {a_sqrt}")
print(f"Rydberg exp(a): {a} eV")
print(f"b_sqrt: {b_sqrt}")
print(f"b: {b}")

x = np.linspace(24,36,100)

plt.plot(x, linear(x,params[0],params[1]), label="Linearer Fit")
plt.plot(Z, energie_sqrt, "x", color = "red", label = "Messwerte")
plt.xlabel(r"$z_{\text{eff}}$")
plt.ylabel(r"$\sqrt{E_{\text{K}}}\,\left[ \sqrt{\text{eV}}\right]$")
plt.xlim([24,36])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Rydberg.pdf")
plt.close()

# Diskussion
def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100)

sigma_1_abweichung = rel_Abweichung (sigma_1_exp, sigma_1_theo)
sigma_2_abweichung = rel_Abweichung (sigma_2_exp, sigma_2_theo)
sigma_3_abweichung = rel_Abweichung (sigma_3_exp, sigma_3_theo)

print("----------------")
print("rel. Abweichungen")
print("----------------")
print(f"sigma_1 Abweichung: {sigma_1_abweichung}%")
print(f"sigma_2 Abweichung: {sigma_2_abweichung}%")
print(f"sigma_3 Abweichung: {sigma_3_abweichung}%")

theta_K_Sr_rel = rel_Abweichung(Theta_K_Sr_exp, 11.04)
theta_K_Zr_rel = rel_Abweichung(Theta_K_Zr_exp, 9.86)
theta_K_Br_rel = rel_Abweichung(Theta_K_Br_exp, 13.23)
theta_K_Zn_rel = rel_Abweichung(Theta_K_Zn_exp, 18.60)
theta_K_Ga_rel = rel_Abweichung(Theta_K_Ga_exp, 17.29)

print("----------------")
print(f"Theta_K_Sr Abweichung: {theta_K_Sr_rel}%")
print(f"Theta_K_Zr Abweichung: {theta_K_Zr_rel}%")
print(f"Theta_K_Br Abweichung: {theta_K_Br_rel}%")
print(f"Theta_K_Zn Abweichung: {theta_K_Zn_rel}%")
print(f"Theta_K_Ga Abweichung: {theta_K_Ga_rel}%")

E_Sr_rel = rel_Abweichung(energie_Sr_exp, 16.10*1000)
E_Zr_rel = rel_Abweichung(energie_Zr_exp, 17.99*1000)
E_Br_rel = rel_Abweichung(energie_Br_exp, 13.47*1000)
E_Zn_rel = rel_Abweichung(energie_Zn_exp, 9.65*1000)
E_Ga_rel = rel_Abweichung(energie_Ga_exp, 10.37*1000)

print("----------------")
print(f"E_Sr Abweichung: {E_Sr_rel} %")
print(f"E_Zr Abweichumg: {E_Zr_rel} %")
print(f"E_Br Abweichung: {E_Br_rel} %")
print(f"E_Zn Abweichung: {E_Zn_rel} %")
print(f"E_Ga Abweichung: {E_Ga_rel} %")

sigma_Sr_rel = rel_Abweichung(sigma_Sr_exp, 4.00)
sigma_Zr_rel = rel_Abweichung(sigma_Zr_exp, 4.10)
sigma_Br_rel = rel_Abweichung(sigma_Br_exp, 3.85)
sigma_Zn_rel = rel_Abweichung(sigma_Zn_exp, 3.56)
sigma_Ga_rel = rel_Abweichung(sigma_Ga_exp, 3.61)

print("----------------")
print(f"sigma_Sr Abweichung: {sigma_Sr_rel}%")
print(f"sigma_Zr Abweichung: {sigma_Zr_rel}%")
print(f"sigma_Br Abweichung: {sigma_Br_rel}%")
print(f"sigma_Zn Abweichung: {sigma_Zn_rel}%")
print(f"sigma_Ga Abweichung: {sigma_Ga_rel}%")

R_inf_rel = rel_Abweichung(unp.nominal_values(a), R)
print("----------------")
print(f"R_inf Abweichung: {R_inf_rel}%")
