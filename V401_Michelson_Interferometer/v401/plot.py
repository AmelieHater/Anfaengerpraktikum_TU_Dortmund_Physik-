import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from uncertainties import ufloat

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

lambda_theo = 645 * 10**(-9) # m (Diodenlaser)
u = 5.046
T_0 = 273.15
T = 21.6 + 273.15 # K 
p_0 = 1019.2 * 100 # Pa
# p_0 = 1013.25 * 100 #Pa

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
z1 = ufloat(np.mean(z1_0), np.std(z1_0))
z2 = ufloat(np.mean(z2_0), np.std(z2_0))
print(f'z1 = {z1}')
print(f'z2 = {z2}')

D = 50 * 10**(-3) # m
delta_n1 = (z1 * lambda_theo)/(2*D)
delta_n2 = (z2 * lambda_theo)/(2*D)
print(f'delta_n1 = {delta_n1}')
print(f'delta_n2 = {delta_n2}')

p = 500 * 133.3 #Pa
n_exp1 = 1 + delta_n1 * (T/T_0) * (p_0/p)
n_exp2 = 1 + delta_n2 * (T/T_0) * (p_0/p)
n_theo = 1.000292
print(f'Delta p = {p/100} hPa')
print(f'n_exp1 = {n_exp1}')
print(f'n_exp2 = {n_exp2}')


print('Diskussion')
lambda_rel = rel_Abweichung(lambda_exp, lambda_theo)
n_rel1 = rel_Abweichung(n_exp1, n_theo)
n_rel2 = rel_Abweichung(n_exp2, n_theo)

print(f'rel. Abw. Wellenlänge = {lambda_rel} %')
print(f'rel. Abw. Brechungsindex 1 = {n_rel1}')
print(f'rel. Abw. Brechungsindex 2 = {n_rel2}')
# Umgebungsdruck 1019.2HPa
