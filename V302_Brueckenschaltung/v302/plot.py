import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as sc

# ------------------------ A ------------------------
print("Aufgabe a: Wheatstonsche Brücke")
# Wert 13
R_2_13, R_3_13, R_4_13 = np.genfromtxt("Messwerte/WheatstonscheBrücke_1.txt", unpack=True)

delta_R_34_13 = (R_3_13/R_4_13)* 0.005
R_34_13 = unp.uarray(R_3_13/R_4_13, delta_R_34_13)
R_x_13 = R_2_13 * R_34_13
# print (f"R_x_13: {R_x_13} in ohm")
R_x_13 = np.mean(R_x_13)
print (f"R_x_13: {R_x_13} in ohm")

#Wert 14
R_2_14, R_3_14, R_4_14 = np.genfromtxt("Messwerte/WheatstonscheBrücke_2.txt", unpack=True)

delta_R_34_14 = (R_3_14/R_4_14)* 0.005
R_34_14 = unp.uarray(R_3_14/R_4_14, delta_R_34_14)
R_x_14 = R_2_14 * R_34_14
# print (f"R_x_14: {R_x_14} in ohm")
R_x_14 = np.mean(R_x_14)
print (f"R_x_14: {R_x_14} in ohm")

# ------------------------ B ------------------------
print("Aufgabe b: Kapazitätenmessbrücke")
# Wert 15
C_2_15 ,R_2_15, R_3_15, R_4_15 = np.genfromtxt("Messwerte/Kapazitätsmessbrücke_1.txt", unpack=True)

delta_R_34_15 = (R_3_15/R_4_15)* 0.005
R_34_15 = unp.uarray(R_3_15/R_4_15, delta_R_34_15)
R_x_15 = R_2_15 * R_34_15
# print (f"R_x_15: {R_x_15} in ohm")
R_x_15 = np.mean(R_x_15)
print (f"R_x_15: {R_x_15} in ohm")

delta_R_43_15 = (R_4_15/R_3_15)* 0.005
R_43_15 = unp.uarray(R_4_15/R_3_15, delta_R_43_15)
C_x_15 = C_2_15 * R_43_15
# print (f"C_x_15: {C_x_15} in nF")
C_x_15 = np.mean(C_x_15)
print (f"C_x_15: {C_x_15} in nF")

# Wert 8
C_2_8 ,R_2_8, R_3_8, R_4_8 = np.genfromtxt("Messwerte/Kapazitätsmessbrücke_2.txt", unpack=True)

delta_R_34_8 = (R_3_8/R_4_8)* 0.005
R_34_8 = unp.uarray(R_3_8/R_4_8, delta_R_34_8)
R_x_8 = R_2_8 * R_34_8
# print (f"R_x_8: {R_x_8} in ohm")
R_x_8 = np.mean(R_x_8)
print (f"R_x_8: {R_x_8} in ohm")

delta_R_43_8 = (R_4_8/R_3_8)* 0.005
R_43_8 = unp.uarray(R_4_8/R_3_8, delta_R_43_8)
C_x_8 = C_2_8 * R_43_8
# print (f"C_x_8: {C_x_8} in nF")
C_x_8 = np.mean(C_x_8)
print (f"C_x_8: {C_x_8} in nF")

# ------------------------ C ------------------------
print("Aufgabe c: Induktivitätsbrücke")
# Wert 19
L_2_19, R_2_19, R_3_19, R_4_19 = np.genfromtxt("Messwerte/Induktivitätsbrücke_1.txt", unpack=True)

delta_R_34_19 = (R_3_19/R_4_19)* 0.005
R_34_19 = unp.uarray(R_3_19/R_4_19, delta_R_34_19)

delta_R_43_19 = (R_4_19/R_3_19)* 0.005
R_43_19 = unp.uarray(R_4_19/R_3_19, delta_R_43_19)
# print (f"R_34_19: {R_34_19}")
R_x_19 = R_2_19 * R_34_19
# print (f"R_x_19: {R_x_19} in ohm")
R_x_19 = np.mean(R_x_19)
print (f"R_x_19: {R_x_19} in ohm")

L_x_19 = L_2_19 * R_43_19
# L_x_19 = L_2_19 * R_34_19
# print (f"L_x_19: {L_x_19} in mH")
L_x_19 = np.mean(L_x_19)
print (f"L_x_19: {L_x_19} in mH")

# Wert 16
L_2_16, R_2_16, R_3_16, R_4_16 = np.genfromtxt("Messwerte/Induktivitätsbrücke_2.txt", unpack=True)

delta_R_34_16 = (R_3_16/R_4_16)* 0.005
R_34_16 = unp.uarray(R_3_16/R_4_16, delta_R_34_16)

delta_R_43_16 = (R_4_16/R_3_16)* 0.005
R_43_16 = unp.uarray(R_4_16/R_3_16, delta_R_43_16) 
# R_x_16 = R_2_16 * R_43_16
R_x_16 = R_2_16 * R_34_16
# print (f"R_x_16: {R_x_16} in ohm")
R_x_16 = np.mean(R_x_16)
print (f"R_x_16: {R_x_16} in ohm")

L_x_16 = L_2_16 * R_43_16
# L_x_16 = L_2_16 * R_34_16
# print (f"L_x_16: {L_x_16} in mH")
L_x_16 = np.mean(L_x_16)
print (f"L_x_16: {L_x_16} in mH")

# ------------------------ D ------------------------
print("Aufgabe d: Maxwell-Brücke")
# Wert 19
R_2_19_d, R_3_19_d, R_4_19_d, C_4_19_d = np.genfromtxt("Messwerte/Maxwellbrücke_1.txt", unpack=True)
C_4_19_d = C_4_19_d * 10**(-9)
delta_R_34_19_d = (R_3_19_d/R_4_19_d)* 0.005
R_34_19_d = unp.uarray(R_3_19_d/R_4_19_d, delta_R_34_19_d)
R_x_19_d = R_2_19_d * R_34_19_d
print (f"R_x_19_d: {R_x_19_d} in ohm")
R_x_19_d = np.mean(R_x_19_d)
print (f"R_x_19_d: {R_x_19_d} in ohm")

L_x_19_d = R_2_19_d * R_3_19_d * C_4_19_d
# print (f"L_x_19_d: {L_x_19_d * 10**3} in mH")
L_x_19_d = np.mean(L_x_19_d)* 10**3 # in mH
print (f"L_x_19_d: {L_x_19_d} in mH")

# Wert 16
R_2_16_d, R_3_16_d, R_4_16_d, C_4_16_d = np.genfromtxt("Messwerte/Maxwellbrücke_2.txt", unpack=True)
C_4_16_d = C_4_16_d * 10**(-9)
delta_R_34_16_d = (R_3_16_d/R_4_16_d)* 0.005
R_34_16_d = unp.uarray(R_3_16_d/R_4_16_d, delta_R_34_16_d)
R_x_16_d = R_2_16_d * R_34_16_d
# print (f"R_x_16_d: {R_x_16_d} in ohm")
R_x_16_d = np.mean(R_x_16_d)
print (f"R_x_16_d: {R_x_16_d} in ohm")

L_x_16_d = R_2_16_d * R_3_16_d * C_4_16_d
# print (f"L_x_16_d: {L_x_16_d * 10**3} in mH")
L_x_16_d = np.mean(L_x_16_d)* 10**3 # in mH
print (f"L_x_16_d: {L_x_16_d} in mH")

# ------------------------ E ------------------------
f_Hz, U_mV = np.genfromtxt("Messwerte/WienRobinsonBrücke.txt", unpack=True)

x = f_Hz/400 # normiert auf kleinste Amplitude
y = (U_mV*10**(-3))**2

def Theoriekurve(x):
    return ((x**2 -1)**2 / (9*((1 - x**2)**2 + 9 * x**2)))

x_neu = np.linspace(0.1, 14, 1000)
fig, ax1 = plt.subplots(1, 1)
ax1.plot(x, y, "x", label="Messwerte")
ax1.plot(x_neu, Theoriekurve(x_neu), "-", label="Theoriekurve")
ax1.set_xscale("log")
ax1.set_xlim([0.1,14])
ax1.set_xlabel(r"$\Omega$")
ax1.set_ylabel(r"$\left|\frac{U_{Br}}{U_{S}}\right|^2$")
ax1.legend(loc="best")
ax1.grid()
fig.savefig("plot.pdf")

# ------------------------ F ------------------------
print("Aufgabe f: Klirrfaktor")
f_2_f = np.sqrt(Theoriekurve(2))
U_br_2_f = 15
U_2_f = U_br_2_f/f_2_f
print (f"U_2_f: {U_2_f} in mV")
print (f"f_2_f: {f_2_f}")
k = U_2_f* 10**(-3)
print(f"k: {k}")

# Diskussion -----------------------------------------
def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100)
#a
print("Diskussion:")
print("a:")
R_x_13_theo = 319.5
R_x_14_theo = 900
rel_R_x_13 = rel_Abweichung(R_x_13, R_x_13_theo)
print(f"rel_R_x_13: {rel_R_x_13}, bei exp: {R_x_13} und theo: {R_x_13_theo}")

rel_R_x_14 = rel_Abweichung(R_x_14, R_x_14_theo)
print(f"rel_R_x_14: {rel_R_x_14}, bei exp: {R_x_14} und theo: {R_x_14_theo}")

#b
print("b:")
R_x_15_theo = 473
R_x_8_theo = 564
C_x_15_theo = 652
C_x_8_theo = 294.1
rel_R_x_15 = rel_Abweichung(R_x_15, R_x_15_theo)
print(f"rel_R_x_15: {rel_R_x_15}, bei exp: {R_x_15} und theo: {R_x_15_theo}")
rel_C_x_15 = rel_Abweichung(C_x_15, C_x_15_theo)
print(f"rel_C_x_15: {rel_C_x_15}, bei exp: {C_x_15} und theo: {C_x_15_theo}")

rel_R_x_8 = rel_Abweichung(R_x_8, R_x_8_theo)
print(f"rel_R_x_8: {rel_R_x_8}, bei exp: {R_x_8} und theo: {R_x_8_theo}")
rel_C_x_8 = rel_Abweichung(C_x_8, C_x_8_theo)
print(f"rel_C_x_8: {rel_C_x_8}, bei exp: {C_x_8} und theo: {C_x_8_theo}")

#c
print("c:")
R_x_19_theo = 108.7
R_x_16_theo = 411.2
L_x_19_theo = 26.96
L_x_16_theo = 132.71
rel_R_x_19 = rel_Abweichung(R_x_19, R_x_19_theo)
print(f"rel_R_x_19: {rel_R_x_19}, bei exp: {R_x_19} und theo: {R_x_19_theo}")
rel_L_x_19 = rel_Abweichung(L_x_19, L_x_19_theo)
print(f"rel_L_x_19: {rel_L_x_19}, bei exp: {L_x_19} und theo: {L_x_19_theo}")

rel_R_x_16 = rel_Abweichung(R_x_16, R_x_16_theo)
print(f"rel_R_x_16: {rel_R_x_16}, bei exp: {R_x_16} und theo: {R_x_16_theo}")
rel_L_x_16 = rel_Abweichung(L_x_16, L_x_16_theo)
print(f"rel_L_x_16: {rel_L_x_16}, bei exp: {L_x_16} und theo: {L_x_16_theo}")

#d
print("d:")
rel_R_x_19_d = rel_Abweichung(R_x_19_d, R_x_19_theo)
print(f"rel_R_x_19_d: {rel_R_x_19_d}, bei exp: {R_x_19_d} und theo: {R_x_19_theo}")
rel_L_x_19_d = rel_Abweichung(L_x_19_d, L_x_19_theo)
print(f"rel_L_x_19_d: {rel_L_x_19_d}, bei exp: {L_x_19_d} und theo: {L_x_19_theo}")

rel_R_x_16_d = rel_Abweichung(R_x_16_d, R_x_16_theo)
print(f"rel_R_x_16_d: {rel_R_x_16_d}, bei exp: {R_x_16_d} und theo: {R_x_16_theo}")
rel_L_x_16_d = rel_Abweichung(L_x_16_d, L_x_16_theo)
print(f"rel_L_x_16_d: {rel_L_x_16_d}, bei exp: {L_x_16_d} und theo: {L_x_16_theo}")