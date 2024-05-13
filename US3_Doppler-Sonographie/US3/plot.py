import matplotlib.pyplot as plt
import numpy as np

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent


Pumpgeschwindigkeit_10mm, delta_nu_von_Pumpgeschwindigkeit_10, Winkel_grad_von_Pumpgeschwindigkeit_10 = np.genfromtxt("content/Messwerte/Geschwindigkeit_10_nach_Winkel.txt", unpack = True)
Pumpgeschwindigkeit_16mm, delta_nu_von_Pumpgeschwindigkeit_16, Winkel_grad_von_Pumpgeschwindigkeit_16 = np.genfromtxt("content/Messwerte/Geschwindigkeit_16_nach_Winkel.txt", unpack = True)


#Strömungsgeschwindigkeit als Funktion der Dopplerwinkel für 5 verschiedene Flussgeschwindigkeiten bestimmen 
# 10 mm Rohr
#Dopplerwinkel berechnen 
def Dopplerwinkel(theta):
    alpha = 90 - ((np.arcsin(np.sin((theta/180)*np.pi)* (2/3)))/ np.pi)*180
    return alpha
alpha_15 = Dopplerwinkel(15)
alpha_30 = Dopplerwinkel(30)
alpha_45 = Dopplerwinkel(45)
#print("Dopplerwinkel:",alpha_15, alpha_30, alpha_45)

#Aus Dopplerverschiebung die Strömungsgeschwindigkeit bestimmen 
nu_0 = 2e6 #in Hertz, was ja wahrscheinlich der Unterschied sein wird
c = 1800 #in m/s
def Stroemungsgeschwindigkeit(alpha, nu_0, c, delta_nu):
    v_St = (c * delta_nu)/(2*nu_0*(np.cos((alpha/180)*np.pi)))
    return v_St

delta_nu_von_15_10 = delta_nu_von_Pumpgeschwindigkeit_10[:5]
delta_nu_von_30_10 = delta_nu_von_Pumpgeschwindigkeit_10[5:10]
delta_nu_von_45_10 = delta_nu_von_Pumpgeschwindigkeit_10[10:]


v_St_15_10 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_von_15_10)
v_St_30_10 = Stroemungsgeschwindigkeit(alpha_30, nu_0, c, delta_nu_von_30_10)
v_St_45_10 = Stroemungsgeschwindigkeit(alpha_45, nu_0, c, delta_nu_von_45_10)
#print(v_St_15_10)
#Das gleiche für das 16er Rohr

delta_nu_von_15_16 = delta_nu_von_Pumpgeschwindigkeit_16[:5]
delta_nu_von_30_16 = delta_nu_von_Pumpgeschwindigkeit_16[5:10]
delta_nu_von_45_16 = delta_nu_von_Pumpgeschwindigkeit_16[10:]

v_St_15_16 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_von_15_16)
v_St_30_16 = Stroemungsgeschwindigkeit(alpha_30, nu_0, c, delta_nu_von_30_16)
v_St_45_16 = Stroemungsgeschwindigkeit(alpha_45, nu_0, c, delta_nu_von_45_16)

#Ausgleichsgerade 
x_Ausgleichsgerade = np.concatenate((abs(v_St_15_10), abs(v_St_30_10), abs(v_St_45_10)))
x_Ausgleichgerade_error = np.concatenate((x_Ausgleichsgerade, [0.9]))
y_Ausgleichsgerade = np.concatenate((abs((delta_nu_von_15_10)/(np.cos((alpha_15/180)*np.pi))), abs((delta_nu_von_30_10)/(np.cos((alpha_30/180)*np.pi))), abs((delta_nu_von_45_10)/(np.cos((alpha_45/180)*np.pi)))))
y_Ausgleichgerade_error = np.concatenate((y_Ausgleichsgerade, [1600]))

params, covariance_matrix = np.polyfit(x_Ausgleichsgerade, y_Ausgleichsgerade, deg=1, cov=True)
#print(abs((delta_nu_von_15_10)/(np.cos((alpha_15/180)*np.pi))), abs((delta_nu_von_30_10)/(np.cos((alpha_30/180)*np.pi))), abs((delta_nu_von_45_10)/(np.cos((alpha_45/180)*np.pi))))
errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip("ab", params, errors):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_plot = np.linspace(0.3,0.9, 10)
#Jetzt die ganze Scheiße plotten
fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label="Lineare Regression",
    linewidth=2,
)
ax1.plot(abs(v_St_15_10), abs((delta_nu_von_15_10)/(np.cos((alpha_15/180)*np.pi))), "o", label=r"$\theta = 15°$")
ax1.plot(abs(v_St_30_10), abs((delta_nu_von_30_10)/(np.cos((alpha_30/180)*np.pi))), "o", label=r"$\theta = 30°$")
ax1.plot(abs(v_St_45_10), abs((delta_nu_von_45_10)/(np.cos((alpha_45/180)*np.pi))), "o", label=r"$\theta = 45°$")
ax1.set_xlabel(r"$v \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax1.set_ylabel(r"$\frac{\Delta \nu}{\text{cos}(\alpha)} \, \left[ \text{Hz} \right]$")
ax1.legend(loc="best")

fig.savefig("plot1.pdf")

#print("v_St_15_16: ", v_St_15_16)
#print("v_St_30_16: ", v_St_30_16)
#print("v_St_45_16: ", v_St_45_16)
#Jetzt für das 16 mm Rohr
#Ausgleichsgerade 
y_1 = abs((delta_nu_von_15_16)/(np.cos((alpha_15/180)*np.pi)))
y_2 = abs((delta_nu_von_30_16)/(np.cos((alpha_30/180)*np.pi)))
y_3 = abs((delta_nu_von_45_16)/(np.cos((alpha_45/180)*np.pi)))
#print(y_1, y_2, y_3)
x_Ausgleichsgerade2 = np.concatenate((abs(v_St_15_16), abs(v_St_30_16), abs(v_St_45_16)))
y_Ausgleichsgerade2 = np.concatenate((y_1,y_2, y_3))
params2, covariance_matrix2 = np.polyfit(x_Ausgleichsgerade2, y_Ausgleichsgerade2, deg=1, cov=True)

errors2 = np.sqrt(np.diag(covariance_matrix2))

for name, value, error in zip("cd", params2, errors2):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_plot2 = np.linspace(0.12,0.45, 10)
#Jetzt die ganze Scheiße plotten
fig2, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(
    x_plot2,
    params2[0] * x_plot2 + params2[1],
    label="Lineare Regression",
    linewidth=2,
)
ax1.plot(abs(v_St_15_16), abs((delta_nu_von_15_16)/(np.cos((alpha_15/180)*np.pi))), "o", label=r"$\theta = 15°$")
ax1.plot(abs(v_St_30_16), abs((delta_nu_von_30_16)/(np.cos((alpha_30/180)*np.pi))), "o", label=r"$\theta = 30°$")
ax1.plot(abs(v_St_45_16), abs((delta_nu_von_45_16)/(np.cos((alpha_45/180)*np.pi))), "o", label=r"$\theta = 45°$")
ax1.set_xlabel(r"$v \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax1.set_ylabel(r"$\frac{\Delta \nu}{\text{cos}(\alpha)} \, \left[ \text{Hz} \right]$")
ax1.legend(loc="best")

fig2.savefig("plot2.pdf")
# Apperantly ist die STeigung exakt gleich lol
#a = 11655.012 ± 0.000
#b = 0.000 ± 0.000
#c = 11655.012 ± 0.000
#d = 0.000 ± 0.000

#2. Teil der Auswertung jaja 
Tiefe_mus_10mm, delta_nu_3_10, delta_nu_6_10 = np.genfromtxt("content/Messwerte/Stroemungsprofil_10mm.txt", unpack = True)
Tiefe_mus_16mm, delta_nu_3_16, delta_nu_6_16 = np.genfromtxt("content/Messwerte/Stroemungsprofil_16mm.txt", unpack = True)

Stroemungsgeschwindigkeit_3_10 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_3_10)
Stroemungsgeschwindigkeit_6_10 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_6_10)
Stroemungsgeschwindigkeit_3_16 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_3_16)
Stroemungsgeschwindigkeit_6_16 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_6_16)

print("Stroemungsgeschwindigkeit_3_10: ", Stroemungsgeschwindigkeit_3_10)
print("Stroemungsgeschwindigkeit_6_10: ", Stroemungsgeschwindigkeit_6_10)
print("Stroemungsgeschwindigkeit_3_16: ", Stroemungsgeschwindigkeit_3_16)
print("Stroemungsgeschwindigkeit_6_16: ", Stroemungsgeschwindigkeit_6_16)


fig3, (ax3) = plt.subplots(1, 1, layout="constrained")
ax3.plot(Tiefe_mus_10mm, abs(Stroemungsgeschwindigkeit_3_10), "x", label=r"$\text{Pumpgeschwindigkeit 3 l/min}$")
ax3.plot(Tiefe_mus_10mm, abs(Stroemungsgeschwindigkeit_6_10), "x", label=r"$\text{Pumpgeschwindigkeit 6 l/min}$")
ax3.set_xlabel(r"$\text{Messtiefe}\,\, x\, \left[\mu\text{s}\right]$")
ax3.set_ylabel(r"$\text{Momentangeschwindigkeit}\,\,v \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax3.legend(loc="best")
fig3.savefig("plot3.pdf")

fig4, (ax4) = plt.subplots(1, 1, layout="constrained")
ax4.plot(Tiefe_mus_16mm, abs(Stroemungsgeschwindigkeit_3_16), "x", label=r"$\text{Pumpgeschwindigkeit 3 l/min}$")
ax4.plot(Tiefe_mus_16mm, abs(Stroemungsgeschwindigkeit_6_16), "x", label=r"$\text{Pumpgeschwindigkeit 6 l/min}$")
ax4.set_xlabel(r"$\text{Messtiefe}\,\, x\, \left[\mu\text{s}\right]$")
ax4.set_ylabel(r"$\text{Momentangeschwindigkeit}\,\,v \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax4.legend(loc="best")
fig4.savefig("plot4.pdf")

