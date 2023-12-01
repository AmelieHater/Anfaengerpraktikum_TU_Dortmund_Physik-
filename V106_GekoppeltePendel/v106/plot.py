import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy

T_5_linkes, T_5_rechtes= np.genfromtxt("Laenge1_Einzelnes_Pendel_5T.txt", unpack=True)
T_linkes0 = T_5_linkes/5
T_linkes_mean = np.mean(T_linkes0)
T_linkes_std = scipy.stats.sem(T_linkes0)
T_linkes = ufloat(T_linkes_mean, T_linkes_std)
# print(f"Exp: T_linkes_L1: {T_linkes} in s")

T_rechtes0 = T_5_rechtes/5
T_rechtes_mean = np.mean(T_rechtes0)
T_rechtes_std = scipy.stats.sem(T_rechtes0)
T_rechtes = ufloat(T_rechtes_mean, T_rechtes_std)
# print(f"Exp: T_rechtes_L1: {T_rechtes}in s")

T_5_linkes_L2, T_5_rechtes_L2= np.genfromtxt("Laenge2_Einzelnes_Pendel_5T.txt", unpack=True)
T_linkes0_L2 = T_5_linkes_L2/5
T_linkes_L2_mean = np.mean(T_linkes0_L2)
T_linkes_L2_std = scipy.stats.sem(T_linkes0_L2)
T_linkes_L2 = ufloat(T_linkes_L2_mean, T_linkes_L2_std)
# print(f"Exp: T_linkes_L2: {T_linkes_L2} in s")

T_rechtes0_L2 = T_5_rechtes_L2/5
T_rechtes_L2_mean = np.mean(T_rechtes0_L2)
T_rechtes_L2_std = scipy.stats.sem(T_rechtes0_L2)
T_rechtes_L2 = ufloat(T_rechtes_L2_mean, T_rechtes_L2_std)
# print(f"Exp: T_rechtes_L2: {T_rechtes_L2}in s")

# Länge 1 --------------------------------------------------------------------------
g = 9.81
l_1 = 0.3250 # in m
# Experimetell:
# Länge 1 gleichphasige Schwingung
T_gleich_L1, T_gegen_L1 =  np.genfromtxt("Laenge1_Schwingung.txt", unpack=True)
T_gleich0_L1 = T_gleich_L1/5
T_gleich_L1_mean = np.mean(T_gleich0_L1)
T_gleich_L1_std = scipy.stats.sem(T_gleich0_L1)
T_gleich_L1 = ufloat(T_gleich_L1_mean, T_gleich_L1_std)
# print(f"Exp: T_gleich_L1: {T_gleich_L1} in s")

omega_gleich_L1 = 2 * np.pi * (1/T_gleich_L1)
# print(f"Exp: omega_gleich_L1: {omega_gleich_L1} in 1/s")

# Theoriewerte:
omega_plus_1_theo = unp.sqrt(g/l_1)
print(f"Theo: omega_plus_1_theo: {omega_plus_1_theo} in 1/s")
T_plus_1_theo = 2 * np.pi * unp.sqrt(l_1/g)
print(f"Theo: T_plus_1_theo: {T_plus_1_theo} in s")

# Länge 1 gegenphasige Schwingung
# Experimentell:
T_gegen0_L1 = T_gegen_L1/5
T_gegen_L1_mean = np.mean(T_gegen0_L1)
T_gegen_L1_std = scipy.stats.sem(T_gegen0_L1)
T_gegen_L1 = ufloat(T_gegen_L1_mean, T_gegen_L1_std)
# print(f"Exp: T_gegen_L1: {T_gegen_L1} in s")

omega_gegen_L1 = 2 * np.pi * (1/T_gegen_L1)
# print(f"Exp: omega_gegen_L1: {omega_gegen_L1} in 1/s")

Kopplung_L1 = (T_gleich_L1**2 - T_gegen_L1**2)/(T_gleich_L1**2 + T_gegen_L1**2)
print(f"Kopplungskonstante K_L1: {Kopplung_L1}")

#Theorie:
omega_minus_1_theo = unp.sqrt((g * (1 + Kopplung_L1))/(l_1 * (1 - Kopplung_L1)))
print (f"Theo: omega_minus_1_theo: {omega_minus_1_theo} in 1/s")
T_minus_1_theo = 2 * np.pi * unp.sqrt((l_1 * (1 - Kopplung_L1))/(g * (1 + Kopplung_L1)))
print (f"Theo: T_mius_1_theo: {T_minus_1_theo} in s")

# Länge 1 Gekoppelte Schwingung
# Experimentell:
T_schwingung_L1, T_schwebung_L1 =  np.genfromtxt("Laenge1_Gekoppelte_Schwingung.txt", unpack=True)
T_schwingung0_L1 = T_schwingung_L1/5
T_schwingung_L1_mean = np.mean(T_schwingung0_L1)
T_schwingung_L1_std = scipy.stats.sem(T_schwingung0_L1)
T_schwingung_L1 = ufloat(T_schwingung_L1_mean, T_schwingung_L1_std)
# print(f"Exp: T_schwingung_L1: {T_schwingung_L1} in s")

omega_schwingung_L1 = 2 * np.pi * (1/T_schwingung_L1) #how?????
print(f"Exp?????: omega_schwingung_L1: {omega_schwingung_L1} in 1/s") 

# Theorie: wenn man einen theoretischen Wert ausrechnen kann




# Experimentell:
T_schwebung0_L1 = T_schwebung_L1/5
T_schwebung_L1_mean = np.mean(T_schwebung0_L1)
T_schwebung_L1_std = scipy.stats.sem(T_schwebung0_L1)
T_schwebung_L1 = ufloat(T_schwebung_L1_mean, T_schwebung_L1_std)
# print(f"Exp: T_schwebung_L1: {T_schwebung_L1} in s")

omega_schwebung_L1 = omega_gleich_L1 - omega_gegen_L1
omega_schwebung_test = 2* np.pi * (1/T_schwebung_L1)
print(f"Exp: omega_schwebung_L1: {omega_schwebung_L1} in 1/s")
print(f"Exp: omega_schwebung_L1 test: {omega_schwebung_test} in 1/s")
# Theorie:
T_schwebung_L1_theo = (T_plus_1_theo * T_minus_1_theo) / (T_plus_1_theo - T_minus_1_theo)
print(f"Theo: T_schwebung_L1_theo: {T_schwebung_L1_theo} in s")
omega_schwebung_L1_theo = omega_plus_1_theo - omega_minus_1_theo
print(f"Theo: omega_schwebung_L1_theo: {omega_schwebung_L1_theo} in 1/s")

# Länge 2 --------------------------------------------------------------------------
l_2 = 0.6530 # in m
T_gleich_L2, T_gegen_L2 =  np.genfromtxt("Laenge2_Schwingung.txt", unpack=True)
# Länge 2 gleichphasige Schwingung
# Experimentell:
T_gleich0_L2 = T_gleich_L2/5
T_gleich_L2_mean = np.mean(T_gleich0_L2)
T_gleich_L2_std = scipy.stats.sem(T_gleich0_L2)
T_gleich_L2 = ufloat(T_gleich_L2_mean, T_gleich_L2_std)
# print(f"Exp: T_gleich_L2: {T_gleich_L2} in s")

omega_gleich_L2 = 2 * np.pi * (1/T_gleich_L2)
# print(f"Exp: omega_gleich_L2: {omega_gleich_L2} in 1/s")

# Theorie:
omega_plus_2_theo = unp.sqrt(g/l_2)
print(f"Theo: omega_plus_2_theo: {omega_plus_2_theo} in 1/s")
T_plus_2_theo = 2 * np.pi * unp.sqrt(l_2/g)
print(f"Theo: T_plus_2_theo: {T_plus_2_theo} in s")

# Länge 2 gegenphasige Schwingung
# Experimentell:
T_gegen0_L2 = T_gegen_L2/5
T_gegen_L2_mean = np.mean(T_gegen0_L2)
T_gegen_L2_std = scipy.stats.sem(T_gegen0_L2)
T_gegen_L2 = ufloat(T_gegen_L2_mean, T_gegen_L2_std)
# print(f"Exp: T_gegen_L2: {T_gegen_L2} in s")

omega_gegen_L2 = 2 * np.pi * (1/T_gegen_L2)
# print(f"Exp: omega_gegen_L2: {omega_gegen_L2} in 1/s")

Kopplung_L2 = (T_gleich_L2**2 - T_gegen_L2**2)/(T_gleich_L2**2 + T_gegen_L2**2)
print(f"Kopplungskonstante K_L2: {Kopplung_L2}")

#Theorie:
omega_minus_2_theo = unp.sqrt((g * (1 + Kopplung_L2))/(l_2 * (1 - Kopplung_L2)))
print (f"Theo: omega_minus_2_theo: {omega_minus_2_theo} in 1/s")
T_minus_2_theo = 2 * np.pi * unp.sqrt((l_2 * (1 - Kopplung_L2))/(g * (1 + Kopplung_L2)))
print (f"Theo: T_mius_2_theo: {T_minus_2_theo} in s")

# Länge 2 Gekoppelte Schwingung
# Experimentell:
T_schwingung_L2, T_schwebung_L2 =  np.genfromtxt("Laenge2_Gekoppelte_Schwingung.txt", unpack=True)
T_schwingung0_L2 = T_schwingung_L2/5
T_schwingung_L2_mean = np.mean(T_schwingung0_L2)
T_schwingung_L2_std = scipy.stats.sem(T_schwingung0_L2)
T_schwingung_L2 = ufloat(T_schwingung_L2_mean, T_schwingung_L2_std)
# print(f"Exp: T_schwingung_L2: {T_schwingung_L2} in s")

omega_schwingung_L2 = 2 * np.pi * (1/T_schwingung_L2) #how?????
print(f"Exp?????: omega_schwingung_L2: {omega_schwingung_L2} in 1/s") 

# Theorie: wenn man einen theoretischen Wert ausrechnen kann




# Experimentell:
T_schwebung0_L2 = T_schwebung_L2/5
T_schwebung_L2_mean = np.mean(T_schwebung0_L2)
T_schwebung_L2_std = scipy.stats.sem(T_schwebung0_L2)
T_schwebung_L2 = ufloat(T_schwebung_L2_mean, T_schwebung_L2_std)
# print(f"Exp: T_schwebung_L2: {T_schwebung_L2} in s")

omega_schwebung_L2 = omega_gleich_L2 - omega_gegen_L2
omega_schwebung_test2 = 2 * np.pi * (1/T_schwebung_L2)
print(f"Exp: omega_schwebung_L2: {omega_schwebung_L2} in 1/s")
print(f"Exp: omega_schwebung_L2 test: {omega_schwebung_test2} in 1/s")
# Theorie:
T_schwebung_L2_theo = (T_plus_2_theo * T_minus_2_theo) / (T_plus_2_theo - T_minus_2_theo)
print(f"Theo: T_schwebung_L2_theo: {T_schwebung_L2_theo} in s")
omega_schwebung_L2_theo = omega_plus_2_theo - omega_minus_2_theo
print(f"Theo: omega_schwebung_L2_theo: {omega_schwebung_L2_theo} in 1/s")

#  Abweichungen bzw. rel. Messfehler
# L1:
print("Diskussion \nRelative Fehler L1:")
rel_T_gleich_L1 = abs(T_gleich_L1 - T_plus_1_theo)/ T_plus_1_theo
rel_omega_gleich_L1 = abs(omega_gleich_L1 - omega_plus_1_theo) / omega_plus_1_theo
print (f"rel. T+ L1: {rel_T_gleich_L1}")
print(f"rel. omega+ L1: {rel_omega_gleich_L1}")
rel_T_gegen_L1 = abs(T_gegen_L1 - T_minus_1_theo) / T_minus_1_theo
rel_omega_gegen_L1 = abs(omega_gegen_L1 - omega_minus_1_theo) / omega_minus_1_theo
print (f"rel. T- L1: {rel_T_gegen_L1}")
print(f"rel. omega- L1: {rel_omega_gegen_L1}")
rel_T_schwebung_L1 = abs(T_schwebung_L1 - T_schwebung_L1_theo) / T_schwebung_L1_theo
rel_omega_schwebung_L1 = abs(omega_schwebung_L1 - omega_schwebung_L1_theo) / omega_schwebung_L1_theo
rel_omega_schwebung_L1test = abs(omega_schwebung_test - omega_schwebung_L1_theo) / omega_schwebung_L1_theo
rel_omega_schwebung_L1test2 = abs(omega_schwebung_test - abs(omega_schwebung_L1_theo)) / abs(omega_schwebung_L1_theo)
print (f"rel. T_S L1: {rel_T_schwebung_L1}")
print(f"rel. omega_S L1: {rel_omega_schwebung_L1}")
print(f"rel. omega_S L1 test: {rel_omega_schwebung_L1test}")
print(f"rel. omega_S L1 test2: {rel_omega_schwebung_L1test2}")

#L2
print("Diskussion \nRelative Fehler L2:")
rel_T_gleich_L2 = abs(T_gleich_L2 - T_plus_2_theo)/ T_plus_2_theo
rel_omega_gleich_L2 = abs(omega_gleich_L2 - omega_plus_2_theo) / omega_plus_2_theo
print (f"rel. T+ L2: {rel_T_gleich_L2}")
print(f"rel. omega+ L2: {rel_omega_gleich_L2}")
rel_T_gegen_L2 = abs(T_gegen_L2 - T_minus_2_theo) / T_minus_2_theo
rel_omega_gegen_L2 = abs(omega_gegen_L2 - omega_minus_2_theo) / omega_minus_2_theo
print (f"rel. T- L2: {rel_T_gegen_L2}")
print(f"rel. omega- L2: {rel_omega_gegen_L2}")
rel_T_schwebung_L2 = abs(T_schwebung_L2 - T_schwebung_L2_theo) / T_schwebung_L2_theo
rel_omega_schwebung_L2 = abs(omega_schwebung_L2 - omega_schwebung_L2_theo) / omega_schwebung_L2_theo
rel_omega_schwebung_L2test = abs(omega_schwebung_test2 - omega_schwebung_L2_theo) / omega_schwebung_L2_theo
rel_omega_schwebung_L2test2 = abs(omega_schwebung_test2 - abs(omega_schwebung_L2_theo)) / abs(omega_schwebung_L2_theo)
print (f"rel. T_S L2: {rel_T_schwebung_L2}")
print(f"rel. omega_S L2: {rel_omega_schwebung_L2}")
print(f"rel. omega_S L2 test: {rel_omega_schwebung_L2test}")
print(f"rel. omega_S L2 test2: {rel_omega_schwebung_L2test2}")