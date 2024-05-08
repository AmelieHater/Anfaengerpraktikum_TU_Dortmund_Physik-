import matplotlib.pyplot as plt
import numpy as np

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent


Pumpgeschwindigkeit_10mm, delta_nu_von_Pumpgeschwindigkeit_10, Winkel_grad_von_Pumpgeschwindigkeit_10 = np.genfromtxt("content/Messwerte/Geschwindigkeit_10_nach_Winkel.txt", unpack = True)
Pumpgeschwindigkeit_16mm, delta_nu_von_Pumpgeschwindigkeit_16, Winkel_grad_von_Pumpgeschwindigkeit_16 = np.genfromtxt("content/Messwerte/Geschwindigkeit_16_nach_Winkel.txt", unpack = True)

Tiefe_mus_10mm, delta_nu_3_10, delta_nu_6_10 = np.genfromtxt("content/Messwerte/Stroemungsprofil_10mm.txt", unpack = True)
Tiefe_mus_16mm, delta_nu_3_16, delta_nu_6_16 = np.genfromtxt("content/Messwerte/Stroemungsprofil_16mm.txt", unpack = True)

#Strömungsgeschwindigkeit als Funktion der Dopplerwinkel für 5 verschiedene Flussgeschwindigkeiten bestimmen 
# 10 mm Rohr
#Dopplerwinkel berechnen 
def Dopplerwinkel(theta):
    alpha = 90 - ((np.arcsin(np.sin((theta/180)*np.pi)* (2/3)))/ np.pi)*180
    return alpha
alpha_15 = Dopplerwinkel(15)
alpha_30 = Dopplerwinkel(30)
alpha_45 = Dopplerwinkel(45)
#print(alpha_15, alpha_30, alpha_45)

#Aus Dopplerverschiebung die Strömungsgeschwindigkeit bestimmen 
nu_0 = 2e06
c = 343.2 #in m/s
def Stroemungsgeschwindigkeit(alpha, nu_0, c, delta_nu):
    v_St = (c * delta_nu)/(2*nu_0*(np.cos(alpha)/180)*np.pi)
    return v_St

delta_nu_von_15_10 = delta_nu_von_Pumpgeschwindigkeit_10[:5]
delta_nu_von_30_10 = delta_nu_von_Pumpgeschwindigkeit_10[5:10]
delta_nu_von_45_10 = delta_nu_von_Pumpgeschwindigkeit_10[10:]


v_St_15_10 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_von_15_10)
v_St_30_10 = Stroemungsgeschwindigkeit(alpha_30, nu_0, c, delta_nu_von_30_10)
v_St_45_10 = Stroemungsgeschwindigkeit(alpha_45, nu_0, c, delta_nu_von_45_10)
print(v_St_15_10)
#Das gleiche für das 16er Rohr

delta_nu_von_15_16 = delta_nu_von_Pumpgeschwindigkeit_16[:5]
delta_nu_von_30_16 = delta_nu_von_Pumpgeschwindigkeit_16[5:10]
delta_nu_von_45_16 = delta_nu_von_Pumpgeschwindigkeit_16[10:]

v_St_15_16 = Stroemungsgeschwindigkeit(alpha_15, nu_0, c, delta_nu_von_15_16)
v_St_30_16 = Stroemungsgeschwindigkeit(alpha_30, nu_0, c, delta_nu_von_30_16)
v_St_45_16 = Stroemungsgeschwindigkeit(alpha_45, nu_0, c, delta_nu_von_45_16)

#Ausgleichsgerade 
x_Ausgleichsgerade = np.concatenate((abs(v_St_15_10), abs(v_St_30_10), abs(v_St_45_10)))
y_Ausgleichsgerade = np.concatenate((abs((delta_nu_von_15_10)/((np.cos(alpha_15)/180)*np.pi)), abs((delta_nu_von_30_10)/((np.cos(alpha_30)/180)*np.pi)), abs((delta_nu_von_45_10)/((np.cos(alpha_45)/180)*np.pi))))
params, covariance_matrix = np.polyfit(x_Ausgleichsgerade, y_Ausgleichsgerade, deg=1, cov=True)

errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip("ab", params, errors):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_plot = np.linspace(3, 30.5, 10)
#Jetzt die ganze Scheiße plotten
fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label="Lineare Regression",
    linewidth=2,
)
ax1.plot(abs(v_St_15_10), abs((delta_nu_von_15_10)/((np.cos(alpha_15)/180)*np.pi)), "o", label=r"$\theta = 15°$")
ax1.plot(abs(v_St_30_10), abs((delta_nu_von_30_10)/((np.cos(alpha_30)/180)*np.pi)), "o", label=r"$\theta = 30°$")
ax1.plot(abs(v_St_45_10), abs((delta_nu_von_45_10)/((np.cos(alpha_45)/180)*np.pi)), "o", label=r"$\theta = 45°$")
ax1.set_xlabel(r"$v_{\text{St}} \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax1.set_ylabel(r"$\frac{\Delta \nu}{\text{cos}(\alpha)} \, \left[ \text{Hz} \right]$")
ax1.legend(loc="best")

fig.savefig("plot1.pdf")

#Jetzt für das 16 mm Rohr
#Ausgleichsgerade 
x_Ausgleichsgerade2 = np.concatenate((abs(v_St_15_16), abs(v_St_30_16), abs(v_St_45_16)))
y_Ausgleichsgerade2 = np.concatenate((abs((delta_nu_von_15_16)/((np.cos(alpha_15)/180)*np.pi)), abs((delta_nu_von_30_16)/((np.cos(alpha_30)/180)*np.pi)), abs((delta_nu_von_45_16)/((np.cos(alpha_45)/180)*np.pi))))
params2, covariance_matrix2 = np.polyfit(x_Ausgleichsgerade2, y_Ausgleichsgerade2, deg=1, cov=True)

errors2 = np.sqrt(np.diag(covariance_matrix2))

for name, value, error in zip("cd", params2, errors2):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_plot2 = np.linspace(1, 16, 10)
#Jetzt die ganze Scheiße plotten
fig2, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(
    x_plot2,
    params2[0] * x_plot2 + params2[1],
    label="Lineare Regression",
    linewidth=2,
)
ax1.plot(abs(v_St_15_16), abs((delta_nu_von_15_16)/((np.cos(alpha_15)/180)*np.pi)), "o", label=r"$\theta = 15°$")
ax1.plot(abs(v_St_30_16), abs((delta_nu_von_30_16)/((np.cos(alpha_30)/180)*np.pi)), "o", label=r"$\theta = 30°$")
ax1.plot(abs(v_St_45_16), abs((delta_nu_von_45_16)/((np.cos(alpha_45)/180)*np.pi)), "o", label=r"$\theta = 45°$")
ax1.set_xlabel(r"$v_{\text{St}} \,  \left[ \frac{\text{m}}{\text{s}} \right]$")
ax1.set_ylabel(r"$\frac{\Delta \nu}{\text{cos}(\alpha)} \, \left[ \text{Hz} \right]$")
ax1.legend(loc="best")

fig2.savefig("plot2.pdf")
# 
#