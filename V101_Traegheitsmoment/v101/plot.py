import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy

#Bestimmung Apperatekonstante D
Federwaagen_Abstand = 0.19975 #in m
phi_D, F_D = np.genfromtxt("D_Brchnng.txt", unpack=True)
#print(F_D * Federwaagen_Abstand)
D_0 = (F_D * Federwaagen_Abstand)/phi_D #in Nm
#print(D_0)
D_mean = np.mean(D_0)
#print(D_mean)
D_std = scipy.stats.sem(D_0)
D = ufloat(D_mean, D_std)
#print("D ist: ",D)

#Eigentraegheitsmoment der Drillachse 
Abstand_a, T5_Periodendauer= np.genfromtxt("I_D_Brchnng.txt", unpack=True)
#print(Abstand_a)
a_quadrat = np.square(Abstand_a)
T_quadrat = np.square(T5_Periodendauer/5)
#print(a_quadrat)
#print(T_quadrat)


x = a_quadrat
y = T_quadrat

x_plot= np.linspace(-0.01, 0.1)
params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
#errors = np.sqrt(np.diag(covariance_matrix))

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(x, y, "o", label="Messwerte")
ax.plot(x_plot, params[0] * x_plot + params[1],
       label="Lineare Regression")
ax.set_xlabel(r"Abstand$\,\,a^2 \,\,[$m$^2]$")
ax.set_ylabel(r"Periodendauer $\,\,T^2\,\,[$s$^2]$")
ax.legend(loc="best")
fig.savefig("plot.pdf")
#print("m ist: ", params[0], " n ist: ", params[1])
#Lineare regression: y = ax+b a = params[0] b = params[1]
#Durchmesser der Gewichte berechnen mit Fehler
Blaues_Hoehe0, Blaues_Durchmesser0, Rotes_Hoehe0, Rotes_Durchmesser0 = np.genfromtxt("Gewichte_Maße.txt", unpack=True)
Gewichte_Durchmesser0 = (Blaues_Durchmesser0 + Rotes_Durchmesser0)/2
Gewichte_Durchmesser_mean = np.mean(Gewichte_Durchmesser0)

Gewichte_Durchmesser_std = scipy.stats.sem(Gewichte_Durchmesser0)
Gewichte_Durchmesser = ufloat(Gewichte_Durchmesser_mean, Gewichte_Durchmesser_std)
#print("Durchmesser der Gewichte ist ",Gewichte_Durchmesser)
#print("Radius der Gewichte ist ",Gewichte_Durchmesser/2)
#Höhe der Gewichte ausrechnen mit Fehler 
Gewicht_Hoehe0 = (Blaues_Hoehe0 + Rotes_Hoehe0) /2
Gewicht_Hoehe_mean = np.mean(Gewicht_Hoehe0)
Gewicht_Hoehe_std = scipy.stats.sem(Gewicht_Hoehe0)
Gewicht_Hoehe = ufloat(Gewicht_Hoehe_mean, Gewicht_Hoehe_std)
#print("Hoehe der Gewichte ist: ",Gewicht_Hoehe)
#print("b ist: ", params[1])

I_D = (params[1]*D) / (4 * np.pi**2) - 2 * 0.26115 * (np.square(Gewichte_Durchmesser/2)/4 + (1/12) * np.square(Gewicht_Hoehe))
#print("I_D ist ",I_D)

#Trägkeitsmoment der Kugel ausrechnen 
#Literaturwert
m_Kugel = 1.1742
Kugel_Durchmesser0 = np.genfromtxt("Kugel.txt", unpack=True)
Kugel_Durchmesser_mean = np.mean(Kugel_Durchmesser0)
Kugel_Durchmesser_std = scipy.stats.sem(Kugel_Durchmesser0)
Kugel_Durchmesser = ufloat(Kugel_Durchmesser_mean, Kugel_Durchmesser_std)
#print("Kugeldurchmesser: ", Kugel_Durchmesser )
I_K_Literatur = (2/5) * m_Kugel * (Kugel_Durchmesser/2)**2 
#print("I_K Literaturwert: ", I_K_Literatur)

#gemessener Wert
T5_Kugel0, T5_Scheibe0 = np.genfromtxt("Periodendauer_Scheibe_Kugel.txt", unpack=True)
T_Kugel0 = T5_Kugel0/5
T_Kugel_mean = np.mean(T_Kugel0)
T_Kugel_std = scipy.stats.sem(T_Kugel0)
T_Kugel = ufloat(T_Kugel_mean, T_Kugel_std)
I_Kugel_gemessen = (np.square(T_Kugel) * D)/(4 * np.pi**2)
#print("T der Kugel: ",T_Kugel)
I_K_gemessen = I_Kugel_gemessen #- I_D #ohne das, da das Drillding schon in D drin ist???
print("I_K gemessen (ohne I_D): ", I_K_gemessen)

#Literaturwert der Scheibe ausrechnen
m_Scheibe = 0.4237 
Scheibe_Hoehe0, Scheibe_Durchmesser0 = np.genfromtxt("Scheibe.txt", unpack=True)
Scheibe_Durchmesser_mean = np.mean(Scheibe_Durchmesser0)
Scheibe_Durchmesser_std = scipy.stats.sem(Scheibe_Durchmesser0)
Scheibe_Durchmesser = ufloat(Scheibe_Durchmesser_mean, Scheibe_Durchmesser_std)
I_S_Literatur = (1/2) * m_Scheibe * (Scheibe_Durchmesser/2)**2
#print("Literaturwert der Scheibe" , I_S_Literatur)

#Trägkeitsmoment der Scheibe ausrechnen 
T_Scheibe0 = T5_Scheibe0/5
T_Scheibe_mean = np.mean(T_Scheibe0)
T_Scheibe_std = scipy.stats.sem(T_Scheibe0)
T_Scheibe = ufloat(T_Scheibe_mean, T_Scheibe_std)
#print("T der Scheibe: ",T_Scheibe)
I_Scheibe_gemessen = (np.square(T_Scheibe) * D)/(4 * np.pi**2)
I_S_gemessen = I_Scheibe_gemessen #- I_D #ohne das, da das Drillding schon in D drin ist???
print("I_S gemessen (ohne I_D): ", I_S_gemessen)

#Trägheitsmoment der 1. Puppenstellung berechnen 
Kopf_Hoehe0 , Arme_Hoehe0 , Torso_Hoehe0 , Beine_Hoehe0 = np.genfromtxt("Puppe_Höhe.txt", unpack=True)
Kopf_Durchmesser00 , Arme_Durchmesser0 , Torso_Durchmesser0 , Beine_Durchmesser0  = np.genfromtxt("Puppe_Durchmesser.txt", unpack=True)
Kopf_Durchmesser0 = Kopf_Durchmesser00[0:5]
#print(Kopf_Durchmesser)

# Mittelung der einzelnen Messerte
Kopf_Hoehe_mean = np.mean(Kopf_Hoehe0)
Kopf_Hoehe_std = scipy.stats.sem(Kopf_Hoehe0)
Kopf_Hoehe = ufloat(Kopf_Hoehe_mean, Kopf_Hoehe_std)
Kopf_Durchmesser_mean = np.mean(Kopf_Durchmesser0)
Kopf_Durchmesser_std = scipy.stats.sem(Kopf_Durchmesser0)
Kopf_Durchmesser = ufloat(Kopf_Durchmesser_mean, Kopf_Durchmesser_std)
Kopf_Durchmesser = ufloat(Kopf_Durchmesser_mean, Kopf_Durchmesser_std)

Arme_Hoehe_mean = np.mean(Arme_Hoehe0)
Arme_Hoehe_std = scipy.stats.sem(Arme_Hoehe0)
Arme_Hoehe = ufloat(Arme_Hoehe_mean, Arme_Hoehe_std)
#print("Armhöhe",Arme_Hoehe)
#print("Arme_Durchmesser0: ", Arme_Durchmesser0)
Arme_Durchmesser_mean = np.mean(Arme_Durchmesser0)
Arme_Durchmesser_std = scipy.stats.sem(Arme_Durchmesser0)
Arme_Durchmesser = ufloat(Arme_Durchmesser_mean, Arme_Durchmesser_std)
#print("Arme Durchmesser",Arme_Durchmesser)

Torso_Hoehe_mean = np.mean(Torso_Hoehe0)
Torso_Hoehe_std = scipy.stats.sem(Torso_Hoehe0)
Torso_Hoehe = ufloat(Torso_Hoehe_mean, Torso_Hoehe_std)
Torso_Durchmesser_mean = np.mean(Torso_Durchmesser0)
Torso_Durchmesser_std = scipy.stats.sem(Torso_Durchmesser0)
Torso_Durchmesser = ufloat(Torso_Durchmesser_mean, Torso_Durchmesser_std)

Beine_Hoehe_mean = np.mean(Beine_Hoehe0)
Beine_Hoehe_std = scipy.stats.sem(Beine_Hoehe0)
Beine_Hoehe = ufloat(Beine_Hoehe_mean, Beine_Hoehe_std)
Beine_Durchmesser_mean = np.mean(Beine_Durchmesser0)
Beine_Durchmesser_std = scipy.stats.sem(Beine_Durchmesser0)
Beine_Durchmesser = ufloat(Beine_Durchmesser_mean, Beine_Durchmesser_std)

#Literaturwert der Puppe in 1. Position
#Massen der einzelnen Körperteile berechnen 
Puppe_Gesamtmasse = 0.1696 # in kg
#Annahme, dass Dichte homogen ist im Körper
#Volumenbestimmung der Puppe
Volumen_1Arm = np.pi * (Arme_Durchmesser/2)**2 * Arme_Hoehe
#print("Volumen Arm: ",Volumen_1Arm)
Volumen_Kopf = np.pi * (Kopf_Durchmesser/2)**2 * Kopf_Hoehe
#print("Volumen Kopf: ",Volumen_Kopf)
Volumen_Torso = np.pi * (Torso_Durchmesser/2)**2 * Torso_Hoehe
#print("Volumen Torso: ",Volumen_Torso)
Volumen_1Bein = np.pi * (Beine_Durchmesser/2)**2 * Beine_Hoehe
#print("Volumen Bein: ",Volumen_1Bein)
Volumen_gesamte_Puppe = 2 * Volumen_1Arm + 2 * Volumen_1Bein + Volumen_Kopf + Volumen_Torso
Volumen_gesamte_Puppe_in_cm = Volumen_gesamte_Puppe * 1000000
#print("Volumen der gesamten Puppe in cm^3: ",Volumen_gesamte_Puppe_in_cm)

#Anteil der einzelnen Volumina am Gesamtvolumen
Prozent_Arm = Volumen_1Arm/Volumen_gesamte_Puppe
Masse_Arm = Puppe_Gesamtmasse * Prozent_Arm
Masse_Arme = 2* Masse_Arm
#print("Masse Arm:", Masse_Arm)

Prozent_Bein = Volumen_1Bein/Volumen_gesamte_Puppe
Masse_Bein = Puppe_Gesamtmasse * Prozent_Bein
Masse_Beine = 2* Masse_Bein
#print("Masse Bein:", Masse_Bein)

Prozent_Kopf = Volumen_Kopf/Volumen_gesamte_Puppe
Masse_Kopf = Puppe_Gesamtmasse * Prozent_Kopf
#print("Masse Kopf:", Masse_Kopf)

Prozent_Torso = Volumen_Torso/Volumen_gesamte_Puppe
Masse_Torso = Puppe_Gesamtmasse * Prozent_Torso
#print("Masse Torso:", Masse_Torso)
#print("Gesamtmasse:", Masse_Arme + Masse_Beine + Masse_Kopf + Masse_Torso)

#I_gesamt der Puppe in Position 1 ausrechnen
I_Kopf_Position_1 = (1/2) * Masse_Kopf * (Kopf_Durchmesser/2)**2
I_Torso_Position_1 = (1/2) * Masse_Torso * (Torso_Durchmesser/2)**2
I_Bein_Position_1 = (1/2) * Masse_Bein * (Beine_Durchmesser/2)**2 + Masse_Bein * (Beine_Durchmesser/2)**2
I_Beine_Position_1 = 2 * I_Bein_Position_1
I_Arm_Position_1 = Masse_Arm * ((Arme_Durchmesser**2)/4+(Arme_Hoehe**2)/12) + Masse_Arm * (((Torso_Durchmesser/2) + (Arme_Hoehe)/2))**2
I_Arme_Position_1 = 2 * I_Arm_Position_1
I_gesamt_Position_1 = I_Kopf_Position_1 + I_Torso_Position_1 + I_Beine_Position_1 + I_Arme_Position_1
#print("I_Literatur der Puppe in Position 1: ",I_gesamt_Position_1)

#I_gesamt der Puppe in Position 2 ausrechnen
I_Bein_Position_2 = Masse_Bein * ((Beine_Durchmesser**2)/4+(Beine_Hoehe**2)/12) + Masse_Bein * (Beine_Hoehe/2)**2
I_Beine_Position_2 = 2 * I_Bein_Position_2
I_gesamt_Position_2 = I_Kopf_Position_1 + I_Torso_Position_1 + I_Arme_Position_1 + I_Beine_Position_2
#print("I_Literatur der Puppe in Position 2: ",I_gesamt_Position_2)

#I_gesamt der Puppe in Position 1 durch Messung bestimmen: 
T5_phi_90_Position_1_0, T5_phi_120_Position_1_0 ,T5_phi_90_Position_2_0 ,T5_phi_120_Position_2_0 = np.genfromtxt("Puppe_Periodendauer.txt", unpack=True)

#Mittelwert der Periodendauer von Position 1 bestimmen
T_phi_90_Position_1_0 = T5_phi_90_Position_1_0/5
T_phi_90_Position_1_mean = np.mean(T_phi_90_Position_1_0)
T_phi_90_Position_1_std = scipy.stats.sem(T_phi_90_Position_1_0)
T_phi_90_Position_1 = ufloat(T_phi_90_Position_1_mean, T_phi_90_Position_1_std)
#print("Mittelwert T bei 90 Grad: ", T_phi_90_Position_1)

T_phi_120_Position_1_0 = T5_phi_120_Position_1_0/5
T_phi_120_Position_1_mean = np.mean(T_phi_120_Position_1_0)
T_phi_120_Position_1_std = scipy.stats.sem(T_phi_120_Position_1_0)
T_phi_120_Position_1 = ufloat(T_phi_120_Position_1_mean, T_phi_120_Position_1_std)
#print("Mittelwert T bei 120 Grad: ",T_phi_120_Position_1)

# Da Mittelwert bei 90 und 120 Grad sehr nah beieinander sind scheint T nicht von Phi abzuhängen, weswegen beide Messungen als eine behandelt werden
T_Position_1 = (T_phi_90_Position_1 + T_phi_120_Position_1)/2
#print("Neuer allgemeiner Mittelwert von T ist: ", T_Position_1)
#I wirklich bestimmen
I_Puppe_Position_1_gemessen = (np.square(T_Position_1) * D)/(4 * np.pi**2)
#print("I_ges gemessen: ", I_Puppe_Position_1_gemessen)
#print(I_D)
I_Puppe_Position_1_gemessen_mit_I_D_evtl = I_Puppe_Position_1_gemessen #- I_D #ohne das, da das Drillding schon in D drin ist???
print("I_Puppe ohne I_D abgezoegen Position 1: ",I_Puppe_Position_1_gemessen_mit_I_D_evtl)

#gemessenes Drehmoment der Puppe in 2. Position
#Mittelwert der Periodendauer von Position 2 bestimmen
T_phi_90_Position_2_0 = T5_phi_90_Position_2_0/5
T_phi_90_Position_2_mean = np.mean(T_phi_90_Position_2_0)
T_phi_90_Position_2_std = scipy.stats.sem(T_phi_90_Position_2_0)
T_phi_90_Position_2 = ufloat(T_phi_90_Position_2_mean, T_phi_90_Position_2_std)
#print("Mittelwert T bei 90 Grad Position 2: ", T_phi_90_Position_2)

T_phi_120_Position_2_0 = T5_phi_120_Position_2_0/5
T_phi_120_Position_2_mean = np.mean(T_phi_120_Position_2_0)
T_phi_120_Position_2_std = scipy.stats.sem(T_phi_120_Position_2_0)
T_phi_120_Position_2 = ufloat(T_phi_120_Position_2_mean, T_phi_120_Position_2_std)
#print("Mittelwert T bei 120 Grad Position 2: ",T_phi_120_Position_2)

# Da Mittelwert bei 90 und 120 Grad sehr nah beieinander sind scheint T nicht von Phi abzuhängen, weswegen beide Messungen als eine behandelt werden
T_Position_2 = (T_phi_90_Position_2 + T_phi_120_Position_2)/2
#print("Neuer allgemeiner Mittelwert von T bei Position 2 ist: ", T_Position_2)
#I wirklich bestimmen
I_Puppe_Position_2_gemessen = (np.square(T_Position_2) * D)/(4 * np.pi**2)
#print("I_ges Position 2 gemessen: ", I_Puppe_Position_2_gemessen)
I_Puppe_Position_2_gemessen_mit_I_D_evtl = I_Puppe_Position_2_gemessen #- I_D #ohne das, da das Drillding schon in D drin ist???
print("I_Puppe Position 2 ohne I_D abgezoegen: ",I_Puppe_Position_2_gemessen_mit_I_D_evtl)


