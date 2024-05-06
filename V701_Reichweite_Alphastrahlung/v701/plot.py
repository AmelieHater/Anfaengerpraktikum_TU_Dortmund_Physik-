import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
from uncertainties import ufloat
import uncertainties.unumpy as unp
import random 

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

Druck_1_mbar, Pulszahl_1, max_Anzahl_Channel_1, Channel_max_Anzahl_1 = np.genfromtxt("content/Messwerte/Reichweite.txt", unpack = True)
Druck_2_mbar, Pulszahl_2, max_Anzahl_Channel_2, Channel_max_Anzahl_2 = np.genfromtxt("content/Messwerte/Reichweite_2.txt", unpack = True)
Daten = np.genfromtxt("content/Messwerte/Statistik.txt", unpack = True)

#Berechnung der effektiven Weglänge
x_01 = 4.5e-2 #in meter
x_02 = 6e-2 # in Meter 
p_0 = 1013 # in mbar 
Weglaenge_eff_1 = x_01 * (Druck_1_mbar/p_0)
Weglaenge_eff_2 = x_02 * (Druck_2_mbar/p_0)
#beides in meter

Energie_1_MeV = (4/Channel_max_Anzahl_1[0]) * Channel_max_Anzahl_1
Energie_2_MeV = (4/Channel_max_Anzahl_2[0]) * Channel_max_Anzahl_2

params, covariance_matrix = np.polyfit(Weglaenge_eff_1, Energie_1_MeV, deg=1, cov=True)
errors = np.sqrt(np.diag(covariance_matrix))

params_2, covariance_matrix_2 = np.polyfit(Weglaenge_eff_2, Energie_2_MeV, deg=1, cov=True)
errors_2 = np.sqrt(np.diag(covariance_matrix_2))
for name, value, error in zip("ab", params_2, errors_2):
    print(f"{name} = {value:.3f} ± {error:.3f}")
x_min_1 = Weglaenge_eff_1[0]-np.abs(Weglaenge_eff_1[1]-Weglaenge_eff_1[0])/2
x_max_1 = Weglaenge_eff_1[-1]-np.abs(Weglaenge_eff_1[1]-Weglaenge_eff_1[0])/2

x_min_2 = Weglaenge_eff_2[0]-np.abs(Weglaenge_eff_2[1]-Weglaenge_eff_2[0])/2
x_max_2 = Weglaenge_eff_2[-1]-np.abs(Weglaenge_eff_2[1]-Weglaenge_eff_2[0])/2

x_plot = np.linspace(x_min_1, x_max_1)
x_plot_2 = np.linspace(x_min_2, x_max_2)


fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(Weglaenge_eff_1, Energie_1_MeV, "x", label="Messwerte")
ax1.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label="Lineare Regression",
    linewidth=1,
)
ax1.set_xlabel(r"$\text{effektive Weglänge x} \, \, [\text{m}]$")
ax1.set_ylabel(r"$\text{Energie} \, \,  [\text{MeV}] $")
ax1.set_xlim(x_min_1, x_max_1)
ax1.set_ylim(0.5, 4.5)
ax1.legend(loc="best")
fig.savefig("Plots/plot1.pdf")
#print("Steigung der Reichweite_1: ", params[0])

fig, (ax2) = plt.subplots(1, 1, layout="constrained")
ax2.plot(Weglaenge_eff_2, Energie_2_MeV, "x", label="Messwerte")
ax2.plot(
    x_plot_2,
    params_2[0] * x_plot_2 + params_2[1],
    label="Lineare Regression",
    linewidth=1,
)
ax2.set_xlabel(r"$\text{effektive Weglänge x} \, \, [\text{m}]$")
ax2.set_ylabel(r"$\text{Energie} \, \,  [\text{MeV}] $")
ax2.set_xlim(x_min_2, x_max_2)
ax2.set_ylim(0.5, 4.5)
ax2.legend(loc="best")
fig.savefig("Plots/plot2.pdf")
#print("Steigung der Reichweite_2: ", params_2[0])

linear_approx_x_1 = Weglaenge_eff_1[9:12]  
linear_approx_x_2 = Weglaenge_eff_2[7:10] 

linear_approx_1 = Pulszahl_1[9:12]  
linear_approx_2 = Pulszahl_2[7:10]

x_plot_3 = np.linspace(0.019, 0.026)
x_plot_4 = np.linspace(0.020, 0.027)

params_3, covariance_matrix_3 = np.polyfit(linear_approx_x_1, linear_approx_1, deg=1, cov=True)
errors_3 = np.sqrt(np.diag(covariance_matrix_3))

params_4, covariance_matrix_4 = np.polyfit(linear_approx_x_2, linear_approx_2, deg=1, cov=True)
errors_4 = np.sqrt(np.diag(covariance_matrix_4))

x_punkt_1 = (max(Pulszahl_1)/2 - params_3[1])/params_3[0] 
x_punkt_2 = (max(Pulszahl_2)/2 - params_4[1])/params_4[0] 
print("Schnittpunkt bei der 1. Liste: ", x_punkt_1, "Schnittpunkt bei der 2. Liste: ", x_punkt_2)
print("y-Koordinate der 1.: ", max(Pulszahl_1)/2, "y-Koordinate der 2.: ", max(Pulszahl_2)/2)

x_min_3 = Weglaenge_eff_1[0]-np.abs(Weglaenge_eff_1[1]-Weglaenge_eff_1[0])/2
x_max_3 = Weglaenge_eff_1[-1]+np.abs(Weglaenge_eff_1[1]-Weglaenge_eff_1[0])/2

fig,(ax3) = plt.subplots(1, 1, layout="constrained")
ax3.plot(Weglaenge_eff_1, Pulszahl_1, "x", label="Messwerte")
ax3.plot(
    x_plot_3,
    params_3[0] * x_plot_3 + params_3[1],
    label="Lineare Regression",
    linewidth=1,
)
ax3.hlines(y = max(Pulszahl_1)/2, xmin = x_min_3, xmax = x_max_3, label = "Hälfte der maximalen Pulszahl")
ax3.plot(x_punkt_1, max(Pulszahl_1)/2, "o", color = "green", linewidth=6, label = "Schnittpunkt")
ax3.set_xlabel(r"$\text{effektive Weglänge x} \, \, [\text{m}]$")
ax3.set_ylabel(r"$\text{Anzahl der Pulse} $")
ax3.set_xlim(x_min_3, x_max_3)
ax3.legend(loc="best")
fig.savefig("Plots/plot3.pdf")

x_min_4 = Weglaenge_eff_2[0]-np.abs(Weglaenge_eff_2[1]-Weglaenge_eff_2[0])/2
x_max_4 = Weglaenge_eff_2[-1]+np.abs(Weglaenge_eff_2[1]-Weglaenge_eff_2[0])/2

fig, (ax4) = plt.subplots(1, 1, layout="constrained")
ax4.plot(Weglaenge_eff_2, Pulszahl_2, "x", label="Messwerte")
ax4.plot(
    x_plot_4,
    params_4[0] * x_plot_4 + params_4[1],
    label="Lineare Regression",
    linewidth=1,
)
ax4.hlines(y = max(Pulszahl_2)/2, xmin = x_min_4, xmax = x_max_4, label = "Hälfte der maximalen Pulszahl")
ax4.plot(x_punkt_2, max(Pulszahl_2)/2, "o", color = "green", linewidth=6, label = "Schnittpunkt")
ax4.set_xlabel(r"$\text{effektive Weglänge x} \, \, [\text{m}]$")
ax4.set_ylabel(r"$\text{Anzahl der Pulse} $")
ax4.set_xlim(x_min_4, x_max_4)
ax4.legend(loc="best")
fig.savefig("Plots/plot4.pdf")

#Histogrammversuch 
Fehler_von_Daten = np.sqrt(Daten)
Daten_plus_Fehler = Daten + Fehler_von_Daten
Daten_minus_Fehler = Daten - Fehler_von_Daten

Daten_mit_Fehlern = unp.uarray(Daten, Fehler_von_Daten)
#print(Daten_mit_Fehlern)
Mittelwert = np.mean(Daten_mit_Fehlern)
sigma = np.std(Daten)
print("Mittelwert und Sigma: ", Mittelwert, sigma)
# create a new random number generator with a fixed seed
rng = default_rng(42)
poisson_1 = rng.poisson(Mittelwert.nominal_value, 100) #1. Eingabe ist der Mittelwert 2. ist wie viele Werte 
gauss_1 = rng.normal(Mittelwert.nominal_value, sigma, size = 100)


fig_Gauss, ((ax5, ax6), (ax7, ax8))= plt.subplots(2,2)
ax5.hist(Daten, bins=10, histtype='step', range = [1100,1400], label = 'Daten', linewidth=1.3)
ax5.hist(Daten_plus_Fehler, bins = 10, histtype='step', range = [1100,1400], label = r"$\text{Daten} + \sqrt{N}$", linewidth=1.3, linestyle = "--")
ax5.hist(Daten_minus_Fehler, bins = 10, histtype='step', range = [1060,1400], label = r"$\text{Daten} - \sqrt{N}$", linewidth=1.3, linestyle = "--")
ax5.hist(gauss_1, bins = 10, histtype='step', range = [1060,1400], label = r"$\text{Gaussverteilung}$", linewidth=2)
ax5.set_title("Daten in 10 Bins")
#ax5.legend(loc="best")
ax6.hist(Daten, bins=14, histtype='step', range = [1100,1400], linewidth=1.3)
ax6.hist(Daten_plus_Fehler, bins = 14, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax6.hist(Daten_minus_Fehler, bins = 14, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax6.hist(gauss_1, bins = 14, histtype='step', range = [1060,1400], linewidth=2)
ax6.set_title("Daten in 14 Bins")
#ax6.legend(loc="best")
ax7.hist(Daten, bins=16, histtype='step', range = [1100,1400], linewidth=1.3)
ax7.hist(Daten_plus_Fehler, bins = 16, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax7.hist(Daten_minus_Fehler, bins = 16, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax7.hist(gauss_1, bins = 16, histtype='step', range = [1060,1400], linewidth=2)
ax7.set_title("Daten in 16 Bins")
#ax7.legend(loc="best")
ax8.hist(Daten, bins=18, histtype='step', range = [1100,1400], linewidth=1.3)
ax8.hist(Daten_plus_Fehler, bins = 18, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax8.hist(Daten_minus_Fehler, bins = 18, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax8.hist(gauss_1, bins = 18, histtype='step', range = [1060,1400], linewidth=2)
ax8.set_title("Daten in 18 Bins")
#ax8.legend(loc="best")

fig_Gauss.legend(loc="lower center", ncol=4)
fig_Gauss.tight_layout()
fig_Gauss.subplots_adjust(bottom=0.15)

fig_Gauss.savefig("Plots/plot5.pdf")

#Poisson Verteilung
fig_Poisson, ((ax5, ax6), (ax7, ax8))= plt.subplots(2,2)
ax5.hist(Daten, bins=10, histtype='step', range = [1100,1400], label = 'Daten', linewidth=1.3)
ax5.hist(Daten_plus_Fehler, bins = 10, histtype='step', range = [1100,1400], label = r"$\text{Daten} + \sqrt{N}$", linewidth=1.3, linestyle = "--")
ax5.hist(Daten_minus_Fehler, bins = 10, histtype='step', range = [1060,1400], label = r"$\text{Daten} - \sqrt{N}$", linewidth=1.3, linestyle = "--")
ax5.hist(poisson_1, bins = 10, histtype='step', range = [1060,1400], label = r"$\text{Gaussverteilung}$", linewidth=2)
ax5.set_title("Daten in 10 Bins")
#ax5.legend(loc="best")
ax6.hist(Daten, bins=14, histtype='step', range = [1100,1400], linewidth=1.3)
ax6.hist(Daten_plus_Fehler, bins = 14, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax6.hist(Daten_minus_Fehler, bins = 14, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax6.hist(poisson_1, bins = 14, histtype='step', range = [1060,1400], linewidth=2)
ax6.set_title("Daten in 14 Bins")
#ax6.legend(loc="best")
ax7.hist(Daten, bins=16, histtype='step', range = [1100,1400], linewidth=1.3)
ax7.hist(Daten_plus_Fehler, bins = 16, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax7.hist(Daten_minus_Fehler, bins = 16, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax7.hist(poisson_1, bins = 16, histtype='step', range = [1060,1400], linewidth=2)
ax7.set_title("Daten in 16 Bins")
#ax7.legend(loc="best")
ax8.hist(Daten, bins=18, histtype='step', range = [1100,1400], linewidth=1.3)
ax8.hist(Daten_plus_Fehler, bins = 18, histtype='step', range = [1100,1400], linewidth=1.3, linestyle = "--")
ax8.hist(Daten_minus_Fehler, bins = 18, histtype='step', range = [1060,1400], linewidth=1.3, linestyle = "--")
ax8.hist(poisson_1, bins = 18, histtype='step', range = [1060,1400], linewidth=2)
ax8.set_title("Daten in 18 Bins")
#ax8.legend(loc="best")

fig_Poisson.legend(loc="lower center", ncol=4)
fig_Poisson.tight_layout()
fig_Poisson.subplots_adjust(bottom=0.15)

fig_Poisson.savefig("Plots/plot6.pdf")



print( "Energieverlust 1: ", rel_Abweichung( params[0],  params_2[0]))
print("Energieverlust 2: ",rel_Abweichung( params_2[0],  params[0]))
print("Mittlere Reichweite 1: ", rel_Abweichung( x_punkt_1,  x_punkt_2))
print("Mittlere Reichweite 2: ", rel_Abweichung( x_punkt_2,  x_punkt_1))

#Werte: 
#Steigung der Reichweite_1:  -112.383 ± 2.515
#y-Achsenabschnitt: 3.990 ± 0.039
#Steigung der Reichweite_2:  -117.482 ± 9.825
#y-Achsenabschnitt: 4.189 ± 0.172
#Schnittpunkt bei der 1. Liste:  x = 0.02256093357890863 
#y-Koordinate der 1.: y = 13112.5 
#Schnittpunkt bei der 2. Liste:  x = 0.02302477244638596
#y-Koordinate der 2.: y = 8019.0
#Energieverlust 1:  -4.437097868848706 %
#Energieverlust 2:  -4.643117538183591 %
#Mittlere Reichweite 1:  3.2168918181617276 % 
#Mittlere Reichweite 2:  3.323815362612409 %