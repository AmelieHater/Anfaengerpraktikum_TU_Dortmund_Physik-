import matplotlib.pyplot as plt
import numpy as np

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

# Strom-Spannungs-Kennlinie
U_1, I_1 = np.genfromtxt('content/Messwerte/Strom-Spannungs-Kennlinie.txt', unpack=True)
# Beleuchtungsstärke
U_50, I_50 = np.genfromtxt('content/Messwerte/Beleuchtungsstaerke.txt', unpack=True)


fig1, ax1 = plt.subplots(layout='constrained')
ax1.plot(U_1, I_1, 'x', color='dodgerblue', label='Messdaten')
ax1.plot(U_50, I_50, 'x', color='navy', label='Messdaten bei 50%')
ax1.legend(loc='best')
ax1.set_xlabel(r'$U\,$[V]')
ax1.set_ylabel(r'$I\,$[A]')
ax1.grid()
fig1.savefig('Plots/Plot1.pdf')

fig2, ax2 = plt.subplots(layout='constrained')
ax2.plot(U_1[0:12], np.sqrt(I_1[0:12]), 'x', color='royalblue', label='Messdaten')
ax2.legend(loc='best')
ax2.set_xlabel(r'$U\,$[V]')
ax2.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax2.grid()
fig2.savefig('Plots/Plot2.pdf')

# Spektrale Abhängigkeit
# Grünes Licht
U_G, I_G = np.genfromtxt('content/Messwerte/Gruen.txt', unpack=True)
fig3, ax3 = plt.subplots(layout='constrained')
ax3.plot(U_G, np.sqrt(I_G), 'x', color='green', label='Messdaten')
ax3.legend(loc='best')
ax3.set_xlabel(r'$U\,$[V]')
ax3.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax3.grid()
fig3.savefig('Plots/Plot3.pdf')
# Orangenes Licht
U_O, I_O = np.genfromtxt('content/Messwerte/Orange.txt', unpack=True)
fig4, ax4 = plt.subplots(layout='constrained')
ax4.plot(U_O, np.sqrt(I_O), 'x', color='orangered', label='Messdaten')
ax4.legend(loc='best')
ax4.set_xlabel(r'$U\,$[V]')
ax4.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax4.grid()
fig4.savefig('Plots/Plot4.pdf')
# Violettes Licht
U_V, I_V = np.genfromtxt('content/Messwerte/Violet.txt', unpack=True)
fig5, ax5 = plt.subplots(layout='constrained')
ax5.plot(U_V, np.sqrt(I_V), 'x', color='indigo', label='Messdaten')
ax5.legend(loc='best')
ax5.set_xlabel(r'$U\,$[V]')
ax5.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax5.grid()
fig5.savefig('Plots/Plot5.pdf')
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