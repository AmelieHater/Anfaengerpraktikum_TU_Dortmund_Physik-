import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

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