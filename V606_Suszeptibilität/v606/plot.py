import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as sc
from scipy.stats import sem

x, y = np.genfromtxt("content/Messwerte/Bandpass.txt", unpack=True)


#def gaussian(x, amplitude, mean, stddev):
#   return amplitude * np.exp(-(((x - mean) / 4 / stddev) ** 2))
y_norm = y / 104.7618
#def gaussian(x, a, b):
#   return np.exp(-(x-b)**2*a)
#
x_2 = np.linspace(20, 40, 10000)
#params, covariance_matrix = curve_fit(gaussian, x, y_norm, p0=(31, 13.5))


def Gauss(x, b, a):
    return np.exp(-(x-b)**2*a)
   
params, pcov = curve_fit(Gauss, x, y_norm, p0 = [31.62, 1])            
err = np.sqrt(np.diag(pcov))

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(x, y_norm, "x", label="Kurve")
ax1.plot(x_2, Gauss(x_2, *params))
ax1.hlines(y =(1/np.sqrt(2)), xmin = 31.42, xmax = 31.85, linewidth = 2, linestyle = "-" , color = "green", label = "Güte")
ax1.set_xlabel("f in kHz")
ax1.set_ylabel("U/U_max")
ax1.legend(loc="best")

fig.savefig("plot.pdf")
#print(params)
v_plus = 31.85
v_minus = 31.42
v_0 = 31.61188736
guete = v_0 /(v_plus - v_minus)
print(guete)
# Güte ist 73,51601711627912
# Güte am Gerät war 20 
Guete_Abweichung = abs(20 - 73.51601711627912)/20
print(Guete_Abweichung)
#Abweichung der Güte: 2.675800855813956

#Suszeptibilität berechnen 
# xi = 4 *(F / Q) * (U_Br / U_Sp) mit F: Querschnitt der Spule, 

#Chi Theorie ausrechen
#Für Nd: 0.0029876846900586157
#Für Gd: 0.0136565452217882
#Für Dy: 0.025040857824875584
Chi_Nd = (1.256637062*10**(-6) * (9.27439952 * 10**(-24))**2 * (8/11)**2 * 2.59 * 10**28 * 4.5*5.5)/(3 * 1.380649 * 10**(-23) * 296.15)
#print(Chi_Nd)
Chi_Gd = (1.256637062*10**(-6) * (9.27439952 * 10**(-24))**2 * (2)**2 * 2.46 * 10**28 * 3.5*4.5)/(3 * 1.380649 * 10**(-23) * 296.15)
#print(Chi_Gd)
Chi_Dy = (1.256637062*10**(-6) * (9.27439952 * 10**(-24))**2 * (1.33)**2 * 2.52 * 10**28 * 7.5*8.5)/(3 * 1.380649 * 10**(-23) * 296.15)
#print(Chi_Dy)

#Chi Experiment ausrechnen 
#Q_real ausrechnen 
#Q_real_gd = 1.2119125494921674e-05
#Q_real_nd = 1.6053372250599395e-05
#Q_real_dy = 1.1668289516390781e-05

M_gd = 14.08e-3
M_nd = 18.48e-3
M_dy = 14.38e-3

rho_nd_w = 7.24e3
rho_gd_w = 7.4e3
rho_dy_w = 7.8e3

L_gd = 15.7e-2
L_nd = 15.9e-2
L_dy = 15.8e-2

Q_real_gd = (M_gd)/(L_gd * rho_gd_w)
Q_real_nd = (M_nd)/(L_nd * rho_nd_w)
Q_real_dy = (M_dy)/(L_dy * rho_dy_w)

#print("Q_real:", Q_real_gd)
#print("Q_real:", Q_real_nd)
#print("Q_real:", Q_real_dy)

#Alles mögliche einlesen und Mittelwert bestimmen 
U_Spannung_vorher_Gd, U_Brueckenspannung_nachher_Gd, Widerstand_vorher_minimal_Gd, Widerstand_nachher_minimal_Gd = np.genfromtxt("content/Messwerte/Gd_Werte.txt", unpack = True)
U_Spannung_vorher_Nd, U_Brueckenspannung_nachher_Nd, Widerstand_vorher_minimal_Nd, Widerstand_nachher_minimal_Nd = np.genfromtxt("content/Messwerte/Nd_Werte.txt", unpack = True)
U_Spannung_vorher_Dy, U_Brueckenspannung_nachher_Dy, Widerstand_vorher_minimal_Dy, Widerstand_nachher_minimal_Dy = np.genfromtxt("content/Messwerte/Dy_Werte.txt", unpack = True)

#Spannung vorher
#Spannungen vorher Gd:  31.67+/-0.33
#Spannungen vorher Nd:  34.0+/-1.5
#Spannungen vorher Dy:  34.0+/-0

U_Spannung_vorher_Gd_Mittelwert = ufloat(np.mean(U_Spannung_vorher_Gd), sem(U_Spannung_vorher_Gd))
U_Spannung_vorher_Nd_Mittelwert = ufloat(np.mean(U_Spannung_vorher_Nd), sem(U_Spannung_vorher_Nd))
U_Spannung_vorher_Dy_Mittelwert = ufloat(np.mean(U_Spannung_vorher_Dy), sem(U_Spannung_vorher_Dy))
#print("Spannungen vorher: ",U_Spannung_vorher_Gd_Mittelwert)
#print("Spannungen vorher: ",U_Spannung_vorher_Nd_Mittelwert)
#print("Spannungen vorher: ",U_Spannung_vorher_Dy_Mittelwert)

#Spannung nachher 
#Spannungen nachher Gd:  48.0+/-1.5
#Spannungen nachher Nd:  35.7+/-0.7
#Spannungen nachher Dy:  80.67+/-0.33

U_Brueckenspannung_nachher_Gd_Mittelwert = ufloat(np.mean(U_Brueckenspannung_nachher_Gd), sem(U_Brueckenspannung_nachher_Gd))
U_Brueckenspannung_nachher_Nd_Mittelwert = ufloat(np.mean(U_Brueckenspannung_nachher_Nd), sem(U_Brueckenspannung_nachher_Nd))
U_Brueckenspannung_nachher_Dy_Mittelwert = ufloat(np.mean(U_Brueckenspannung_nachher_Dy), sem(U_Brueckenspannung_nachher_Dy))

#print("Spannungen nachher: ",U_Brueckenspannung_nachher_Gd_Mittelwert)
#print("Spannungen nachher: ",U_Brueckenspannung_nachher_Nd_Mittelwert)
#print("Spannungen nachher: ",U_Brueckenspannung_nachher_Dy_Mittelwert)

# Widerstand vorher minimal
#Widerstand vorher Gd:  554+/-6
#Widerstand vorher Nd:  550+/-5
#Widerstand vorher Dy:  544.7+/-2.8

Widerstand_vorher_minimal_Gd_Mittelwert = ufloat(np.mean(Widerstand_vorher_minimal_Gd), sem(Widerstand_vorher_minimal_Gd))
Widerstand_vorher_minimal_Nd_Mittelwert = ufloat(np.mean(Widerstand_vorher_minimal_Nd), sem(Widerstand_vorher_minimal_Nd))
Widerstand_vorher_minimal_Dy_Mittelwert = ufloat(np.mean(Widerstand_vorher_minimal_Dy), sem(Widerstand_vorher_minimal_Dy))

#print("Widerstand vorher: ",Widerstand_vorher_minimal_Gd_Mittelwert )
#print("Widerstand vorher: ", Widerstand_vorher_minimal_Nd_Mittelwert)
#print("Widerstand vorher: ", Widerstand_vorher_minimal_Dy_Mittelwert)

#Widerstand nachher minimal
#Widerstand nachher Gd:  388+/-10
#Widerstand nachher Nd:  527+/-6
#Widerstand nachher Dy:  225+/-8

Widerstand_nachher_minimal_Gd_Mittelwert = ufloat(np.mean(Widerstand_nachher_minimal_Gd), sem(Widerstand_nachher_minimal_Gd))
Widerstand_nachher_minimal_Nd_Mittelwert = ufloat(np.mean(Widerstand_nachher_minimal_Nd), sem(Widerstand_nachher_minimal_Nd))
Widerstand_nachher_minimal_Dy_Mittelwert = ufloat(np.mean(Widerstand_nachher_minimal_Dy), sem(Widerstand_nachher_minimal_Dy))

#print("Widerstand nachher: ", Widerstand_nachher_minimal_Gd_Mittelwert)
#print("Widerstand nachher: ", Widerstand_nachher_minimal_Nd_Mittelwert)
#print("Widerstand nachher: ", Widerstand_nachher_minimal_Dy_Mittelwert)

#Berechnung nach der 1. Methode 
F = 86.6e-6 #in Quadratmeter

#Chi_experimental Gd:  43.3+/-1.5
#Chi_experimental Nd:  22.6+/-1.1
#Chi_experimental Dy:  70.43+/-0.29

Chi_Gd_experimental_1 = 4 * (F/Q_real_gd) * (U_Brueckenspannung_nachher_Gd_Mittelwert/U_Spannung_vorher_Gd_Mittelwert)
Chi_Nd_experimental_1 = 4 * (F/Q_real_nd) * (U_Brueckenspannung_nachher_Nd_Mittelwert/U_Spannung_vorher_Nd_Mittelwert)
Chi_Dy_experimental_1 = 4 * (F/Q_real_dy) * (U_Brueckenspannung_nachher_Dy_Mittelwert/U_Spannung_vorher_Dy_Mittelwert)

print("Chi_experimental: ", Chi_Gd_experimental_1)
print("Chi_experimental: ", Chi_Nd_experimental_1)
print("Chi_experimental: ", Chi_Dy_experimental_1)

#Berechnung nach der 2. Methode
#Chi_2_experimental Gd:  4.27+/-0.28
#Chi_2_experimental Nd:  0.45+/-0.16
#Chi_2_experimental Dy:  8.70+/-0.23

#Berechnung Delta R (Die Einheiten sind egal, weil in der Formel R/R geteilt wird. Daher ist nur das Verhältnis wichtig)
Delta_R_Gd = abs(Widerstand_nachher_minimal_Gd_Mittelwert - Widerstand_vorher_minimal_Gd_Mittelwert)
Delta_R_Nd = abs(Widerstand_nachher_minimal_Nd_Mittelwert - Widerstand_vorher_minimal_Nd_Mittelwert)
Delta_R_Dy = abs(Widerstand_nachher_minimal_Dy_Mittelwert - Widerstand_vorher_minimal_Dy_Mittelwert)

Chi_Gd_experimental_2 = 2 * (Delta_R_Gd/Widerstand_vorher_minimal_Gd_Mittelwert) * (F/Q_real_gd)
Chi_Nd_experimental_2 = 2 * (Delta_R_Nd/Widerstand_vorher_minimal_Nd_Mittelwert) * (F/Q_real_nd)
Chi_Dy_experimental_2 = 2 * (Delta_R_Dy/Widerstand_vorher_minimal_Dy_Mittelwert) * (F/Q_real_dy)

#print("Chi_2_experimental: ",Chi_Gd_experimental_2)
#print("Chi_2_experimental: ",Chi_Nd_experimental_2)
#print("Chi_2_experimental: ",Chi_Dy_experimental_2)

#Abweichungen berechnen

#Chi 1

Abweichung_chi_1_Gd = Chi_Gd_experimental_1/Chi_Gd
Abweichung_chi_1_Nd = Chi_Nd_experimental_1/Chi_Nd
Abweichung_chi_1_Dy = Chi_Dy_experimental_1/Chi_Dy

print("Abweichung Chi Nr.1: ",Abweichung_chi_1_Gd)
print("Abweichung Chi Nr.1: ",Abweichung_chi_1_Nd)
print("Abweichung Chi Nr.1: ",Abweichung_chi_1_Dy)

#Chi 2

Abweichung_chi_2_Gd = Chi_Gd_experimental_2/Chi_Gd
Abweichung_chi_2_Nd = Chi_Nd_experimental_2/Chi_Nd
Abweichung_chi_2_Dy = Chi_Dy_experimental_2/Chi_Dy

print("Abweichung Nr.2: ",Abweichung_chi_2_Gd)
print("Abweichung Nr.2: ",Abweichung_chi_2_Nd)
print("Abweichung Nr.2: ",Abweichung_chi_2_Dy)

#Abweichung Nr. 1:
#Abweichung Chi Nr.1 Gd:  (3.17+/-0.11)e+03
#Abweichung Chi Nr.1 Nd:  (7.6+/-0.4)e+03
#Abweichung Chi Nr.1 Dy:  2813+/-12

#Abweichung Nr.2:
#Abweichung Nr.2 Gd:  312+/-21
#Abweichung Nr.2 Nd:  (1.5+/-0.5)e+02
#Abweichung Nr.2 Dy:  348+/-9
