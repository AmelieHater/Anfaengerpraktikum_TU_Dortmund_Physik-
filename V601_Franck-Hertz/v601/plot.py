import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

#Mittlere freie Weglänge berechnen: 
T_Celsium = np.array([22.6, 150, 160, 171, 178, 180])
T_Kelvin = T_Celsium + 273.15
w_mittel_cm = (29/55000)*10**(-7)*np.exp((6876/T_Kelvin))
w_mittel_m = w_mittel_cm * 10**(-2)
#print("w_mittel_m ist: ",w_mittel_m) # w_mittel_m ist:  [6.59335831e-03 6.01351227e-06 4.13240884e-06 2.78907059e-06 2.19350942e-06 2.05081157e-06]

#Verhältnis zu a: 
a = 10**(-2)
verhaeltnis_a_w = a / w_mittel_m
#print("verhaeltnis_a_w: ", verhaeltnis_a_w) #verhaeltnis_a_w:  [1.51667777e+00 1.66292169e+03 2.41989609e+03 3.58542377e+03 4.55890452e+03 4.87611839e+03]

#Steigung berechnen: delta_y / delta_x
delta_y_alt, delta_x_alt = np.genfromtxt("content/Messwerte/Steigung_22,6.txt", unpack = True)
#print(delta_x_alt)
#print(delta_y_alt)
delta_y = delta_y_alt * ((62-0.2)/170) #delta_ in nA
delta_x = delta_x_alt * (0.4/10) #delta_x in Volt
#print(delta_x)
#print(delta_y)
Steigung_01 = delta_y/delta_x
U_Kaestchen = np.array([5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175,182.5,187.5,192.5,197.5,202.5,207.5,215,225,235])
U_01 = U_Kaestchen * (0.4/10)

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(U_01, Steigung_01,"x", label="Lokale Steigung")
ax1.set_xlabel(r"$U \left[ \text{V}\right]$")
ax1.set_ylabel(r"$I/U \left[ \text{nA/V}\right]$")
ax1.legend(loc="best")
fig.savefig("plot1.pdf")

#Steigung der 2. Kurve berechnen
delta_y_alt_2, delta_x_alt_2 = np.genfromtxt("content/Messwerte/Steigung_150.txt", unpack = True)
#print(len(delta_x_alt_2))
#print(len(delta_y_alt_2))
delta_y_2 = delta_y_alt_2 * ((5*0.3)/134) #delta_ in nA
delta_x_2 = np.append(delta_x_alt_2[:9] * (2/92), delta_x_alt_2[9:] * (2/88)) #delta_x in Volt
#print(len(delta_y_2))
#print(len(delta_x_2))
#print(delta_x_2)
#print(delta_y_2)
Steigung_01_2 = delta_y_2/delta_x_2
#print(len(Steigung_01_2))
U_Kaestchen_2 = np.array([5,15,25,35,45,55,65,75,86,96,105,115,125,135,145,155,165,175,185,195,205,215,225,235])
U_01_2 = np.append(U_Kaestchen_2[:9] * (2/92), U_Kaestchen_2[9:] * (2/88))
#print(U_01_2)
fig, (ax2) = plt.subplots(1, 1, layout="constrained")
ax2.plot(U_01_2, Steigung_01_2,"x", label="Lokale Steigung")
ax2.set_xlabel(r"$U \left[ \text{V}\right]$")
ax2.set_ylabel(r"$I/U \left[ \text{nA/V}\right]$")
ax2.legend(loc="best")
fig.savefig("plot2.pdf")

#Kontaktpotential berechnen:(durch Kurve bei 22,6 °C)
U_kontakt = 11 - U_01[np.argmax(Steigung_01)] #(U_Beschleunigung - U_höchste_Steigungsänderung)
print(U_kontakt) #2.9000000000000004

#Maxima der Kurve nehmen
U_max_position_kaestchen = np.array([14, 32, 52, 73, 95, 116, 139, 163])
U_max_position_Volt = np.append(U_max_position_kaestchen[:3]*(1/4),U_max_position_kaestchen[3:]*(20/83))
U_dummy = np.linspace(1,8,8)

params, covariance_matrix = np.polyfit(U_dummy, U_max_position_Volt, deg=1, cov=True)

errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip("ab", params, errors):
    print(f"{name} = {value:.3f} ± {error:.3f}")
    #a = 5.096 ± 0.079
    #b = -2.219 ± 0.397

fig, (ax3) = plt.subplots(1, 1, layout="constrained")
ax3.plot(U_dummy, U_max_position_Volt,"x", label="Spannungsposition der Maxima")
ax3.plot(
    U_dummy,
    params[0] * U_dummy + params[1],
    label="Lineare Regression",
    linewidth=3,
)
ax3.set_xlabel(r"$\text{N-tes Maximum}$")
ax3.set_ylabel(r"$ U \left[ \text{V}\right]$")
ax3.legend(loc="best")
fig.savefig("plot3.pdf")

#Berechnung des korresponderienden Wellenlänge
h = 4.1356676969e-15 #eV·s           #6.62607015e-34 #J s
c = 299792458 #m/s
Wellenlaenge = (h * c)/params[0] #params[0] #(h * c)/U * e !
print("Wellenlaenge",Wellenlaenge) # 2.432958053086267e-07








U_diff_max_kaestchen = np.array([18, 20, 21, 23, 21, 23, 24])
U_diff_max_volt = np.append(U_diff_max_kaestchen[:3]*(1/4),U_diff_max_kaestchen[3:]*(20/83))
U_max_kaestchen = np.array([(32+14)/2, (52+32)/2, (73+52)/2, (95+73)/2, (116+95)/2, (139+116)/2, (163+139)/2])
U_max_Volt = np.append(U_max_kaestchen[:3]*(1/4),U_max_kaestchen[3:]*(20/83))
U_max_Volt_mittel = np.mean(U_diff_max_volt)
print("U_max_Volt_mittel: ", U_max_Volt_mittel) # 5.239672977624784
U_max_Volt_mittel_Fehler = sc.stats.sem(U_diff_max_volt)
print(U_max_Volt_mittel_Fehler) #0.16295027898423778


#print(U_max_kaestchen)
#params, covariance_matrix = np.polyfit(U_max_Volt, U_diff_max_volt, deg=1, cov=True)
#
#errors = np.sqrt(np.diag(covariance_matrix))
#
#for name, value, error in zip("ab", params, errors):
#    print(f"{name} = {value:.3f} ± {error:.3f}")
#    #a = 0.033 ± 0.009
#    #b = 4.552 ± 0.216
#
#x_plot = np.linspace(3,38,10)
#fig, (ax3) = plt.subplots(1, 1, layout="constrained")
#ax3.plot(U_max_Volt, U_diff_max_volt,"x", label="Abstand der Maxima zueinander")
#ax3.plot(
#    x_plot,
#    params[0] * x_plot + params[1],
#    label="Lineare Regression",
#    linewidth=3,
#)
#ax3.set_xlabel(r"$U \left[ \text{V}\right]$")
#ax3.set_ylabel(r"$\Delta U \left[ \text{V}\right]$")
#ax3.legend(loc="best")
#fig.savefig("plot3.pdf")

#x = np.linspace(0, 10, 1000)
#y = x ** np.sin(x)
#
#fig, (ax1, ax2) = plt.subplots(1, 2, layout="constrained")
#ax1.plot(x, y, label="Kurve")
#ax1.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
#ax1.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
#ax1.legend(loc="best")
#
#ax2.plot(x, y, label="Kurve")
#ax2.set_xlabel(r"$\alpha \mathbin{/} \unit{\ohm}$")
#ax2.set_ylabel(r"$y \mathbin{/} \unit{\micro\joule}$")
#ax2.legend(loc="best")
#
#fig.savefig("build/plot.pdf")