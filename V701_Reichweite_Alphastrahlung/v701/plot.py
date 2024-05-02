import matplotlib.pyplot as plt
import numpy as np

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

Druck_1_mbar, Pulszahl_1, max_Anzahl_Channel_1, Channel_max_Anzahl_1 = np.genfromtxt("\content\Messwerte\Reichweite.txt", unpack = True)
Druck_2_mbar, Pulszahl_2, max_Anzahl_Channel_2, Channel_max_Anzahl_2 = np.genfromtxt("\content\Messwerte\Reichweite_2.txt", unpack = True)
Daten = np.genfromtxt("\content\Messwerte\Statistik.txt", unpack = True)

#Berechnung der effektiven Weglänge
x_01 = 4.5e-2 #in meter
x_02 = 6e-2 # in Meter 
p_0 = 1013 # in mbar 
Weglaenge_eff_1 = x_01 * (Druck_1_mbar/p_0)
Weglaenge_eff_2 = x_02 * (Druck_2_mbar/p_0)
#beides in meter

Energie_1_mMeV = (4/)

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