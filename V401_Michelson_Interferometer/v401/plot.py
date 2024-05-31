import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from uncertainties import ufloat

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

lambda_theo = 645 * 10**(-9) # m (Diodenlaser)
u = 5.046
T_0 = 21.6 + 273.15 # K 
p_0 = 1019.2 * 10**(-3) # bar

z_0 = np.genfromtxt("content/Messdaten/Messwerte.txt", unpack=True)
z = ufloat(np.mean(z_0), np.std(z_0)) 
# x_1 = ufloat(6, 0.01)
# x_2 =ufloat(11, 0.01)
x = (ufloat(5, 0.01) *10**(-3))/u  # m
lambda_exp = (2 * x) / z # m

print(f'z = {z}')
print(f'x = {x * 10**3} mm')
print(f'Wellenlänge exp = {lambda_exp * 10**9} nm')
print(f'Wellenlänge theo = {lambda_theo * 10**9} nm')

z1_0, z2_0 = np.genfromtxt("content/Messdaten/Vakuum.txt", unpack=True)
z2 = ufloat(np.mean(z2_0), np.std(z2_0))
print(f'z = {z2}')
D = 50 * 10**(-3) # m
delta_n = (z2 * lambda_theo)/(2*D)
print(f'delta_n = {delta_n}')
T = T_0
T_0 = 273.15
p = 20 * 10**(-2) # bar
n_exp = 1 + delta_n * (T/T_0) * (p_0/p)
n_theo = 1.000292
print(f'n_exp = {n_exp}')

print('Diskussion')
lambda_rel = rel_Abweichung(lambda_exp, lambda_theo)
n_rel = rel_Abweichung(n_exp, n_theo)

print(f'rel. Abw. Wellenlänge = {lambda_rel} %')
print(f'rel. Abw. Brechungsindex = {n_rel}')
# Umgebungsdruck 1019.2HPa
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