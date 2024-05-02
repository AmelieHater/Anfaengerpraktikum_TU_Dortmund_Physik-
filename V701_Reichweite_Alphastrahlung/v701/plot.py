import matplotlib.pyplot as plt
import numpy as np

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
#for name, value, error in zip("ab", params, errors):
#    print(f"{name} = {value:.3f} ± {error:.3f}")
x_plot = np.linspace(-0.001, 0.027)
x_plot_2 = np.linspace(-0.001, 0.031)

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(Weglaenge_eff_1, Energie_1_MeV, "x", label="Messwerte")
ax1.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label="Lineare Regression",
    linewidth=1,
)
ax1.set_xlabel(r"$\text{effektive Weglänge x} \, \, [m^{-1}]$")
ax1.set_ylabel(r"$\text{Energie} \, \,  [\text{MeV}] $")
ax1.set_xlim(-0.001, 0.027)
ax1.legend(loc="best")
fig.savefig("plot1.pdf")

fig, (ax2) = plt.subplots(1, 1, layout="constrained")
ax2.plot(Weglaenge_eff_2, Energie_2_MeV, "x", label="Messwerte")
ax2.plot(
    x_plot_2,
    params_2[0] * x_plot_2 + params_2[1],
    label="Lineare Regression",
    linewidth=1,
)
ax2.set_xlabel(r"$\text{effektive Weglänge x} \, \, [m^{-1}]$")
ax2.set_ylabel(r"$\text{Energie} \, \,  [\text{MeV}] $")
ax2.set_xlim(-0.001, 0.031)
ax2.legend(loc="best")
fig.savefig("plot2.pdf")
