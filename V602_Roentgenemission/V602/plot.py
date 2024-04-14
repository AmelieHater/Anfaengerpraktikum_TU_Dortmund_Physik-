import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
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
    return Z - np.sqrt((E/R) - ((alpha**2)*(Z**4)/4))
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


# Strontium
theta2_Sr, r_Sr = np.genfromtxt("content/Messwerte/Strontium.txt", unpack=True)
theta_Sr = theta2_Sr / 2

plt.plot(theta_Sr, r_Sr, "x-", label="Messwerte")
# plt.plot(45.2/2, 4756, "o", color="red", label=r"$K_{\alpha}$")
# plt.plot(40.4/2, 1387, "o", color="orange", label=r"$K_{\beta}$")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([10,12])
plt.ylim([30,120])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Strontium.pdf")
plt.close()

# Zirkonium
theta2_Zr, r_Zr = np.genfromtxt("content/Messwerte/Zirkonium.txt", unpack=True)
theta_Zr = theta2_Zr / 2

plt.plot(theta_Zr, r_Zr, "x-", label="Messwerte")
# plt.plot(45.2/2, 4756, "o", color="red", label=r"$K_{\alpha}$")
# plt.plot(40.4/2, 1387, "o", color="orange", label=r"$K_{\beta}$")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([9,10.9])
plt.ylim([70,180])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Zirkonium.pdf")
plt.close()

# Brom
theta2_Br, r_Br = np.genfromtxt("content/Messwerte/Brom.txt", unpack=True)
theta_Br = theta2_Br / 2

plt.plot(theta_Br, r_Br, "x-", label="Messwerte")
# plt.plot(45.2/2, 4756, "o", color="red", label=r"$K_{\alpha}$")
# plt.plot(40.4/2, 1387, "o", color="orange", label=r"$K_{\beta}$")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([12,13.9])
plt.ylim([10,30])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Brom.pdf")
plt.close()

# Zink
theta2_Zn, r_Zn = np.genfromtxt("content/Messwerte/Zink.txt", unpack=True)
theta_Zn = theta2_Zn / 2

plt.plot(theta_Zn, r_Zn, "x-", label="Messwerte")
# plt.plot(45.2/2, 4756, "o", color="red", label=r"$K_{\alpha}$")
# plt.plot(40.4/2, 1387, "o", color="orange", label=r"$K_{\beta}$")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([19.4,20.9])
plt.ylim([25,430])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Brom.pdf")
plt.close()

# Gallium
theta2_Ga, r_Ga = np.genfromtxt("content/Messwerte/Gallium.txt", unpack=True)
theta_Ga = theta2_Ga / 2

plt.plot(theta_Ga, r_Ga, "x-", label="Messwerte")
# plt.plot(45.2/2, 4756, "o", color="red", label=r"$K_{\alpha}$")
# plt.plot(40.4/2, 1387, "o", color="orange", label=r"$K_{\beta}$")
# plt.axvline(x=45.2/2,color="green" ,linestyle = "dotted")
# plt.axvline(x=40.4/2,color="green" ,linestyle = "dotted")
plt.xlabel(r"$\theta\,$[°]")
plt.ylabel(r"R $\left[\frac{\text{Imp}}{\text{s}}\right]$")
plt.xlim([17,19])
plt.ylim([20,34])
plt.legend(loc="best")
plt.grid()
plt.savefig("content/Plots/Gallium.pdf")
plt.close()

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