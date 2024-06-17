import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy as scipy
import scipy.constants as const
from scipy.optimize import curve_fit

def rel_Abweichung(exp, theo):
    exp = np.abs(unp.nominal_values(exp))
    return (np.abs(exp-theo)/(theo)*100) #ist schon in Prozent

def linear(x, m, b):
    return m*x+b

# Strom-Spannungs-Kennlinie
U_1, I_1 = np.genfromtxt('content/Messwerte/Strom-Spannungs-Kennlinie.txt', unpack=True)
U_1 = -1 * U_1

# Beleuchtungsstärke
U_50, I_50 = np.genfromtxt('content/Messwerte/Beleuchtungsstaerke.txt', unpack=True)
U_50 = -1* U_50

fig1, ax1 = plt.subplots(layout='constrained')
ax1.plot(U_1, I_1, 'x--', color='dodgerblue', label='Messdaten bei ganzer Beleuchtungsstärke')
ax1.plot(U_50, I_50, 'x--', color='navy', label='Messdaten bei halber Beleuchtungsstärke')
ax1.legend(loc='best')
ax1.set_xlabel(r'$U\,$[V]')
ax1.set_ylabel(r'$I\,$[A]')
ax1.set_xlim([-20,2.5])
ax1.grid()
fig1.savefig('Plots/Plot1.pdf')

# Fit für Blaues Licht
U_1_plot = U_1[0:12]
I_1_sqrt_plot = np.sqrt(I_1[0:12])
params1, pcov1 = curve_fit(linear, U_1_plot, I_1_sqrt_plot)
error1 = np.sqrt(np.diag(pcov1))
params1 = unp.uarray(params1, error1)

x1 = np.linspace(0.75, 1.10, 1000)
# x1 = np.linspace(-1.10, -0.75, 1000)
m_blau = unp.nominal_values(params1[0])
b_blau = unp.nominal_values(params1[1])

fig2, ax2 = plt.subplots(layout='constrained')
ax2.plot(U_1_plot, I_1_sqrt_plot, 'x', color='royalblue', label='Messdaten')
ax2.plot(
    x1,
    linear(x1, m_blau, b_blau),
    linestyle = '-',
    color = 'navy',
    label = 'Lineare Regression'
)
ax2.legend(loc='best')
ax2.set_xlabel(r'$U\,$[V]')
ax2.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax2.grid()
ax2.set_ylim([-0.1e-5, 1.55e-5])
ax2.set_xlim([0.75, 1.10])
# ax2.set_xticks(np.arange(-1.10, -0.75, step=0.05))
fig2.savefig('Plots/blau.pdf')

U_G_B = -params1[1]/ params1[0]
print(f'm und b für lin Reg Blau: {params1}')
print(f'Grenspannung Blau: U_G= {U_G_B}V\n mit x=(y-b)/m')

# Spektrale Abhängigkeit
# Grünes Licht
U_G, I_G = np.genfromtxt('content/Messwerte/Gruen.txt', unpack=True)
U_G = -1 * U_G

params2, pcov2 = curve_fit(linear, U_G, np.sqrt(I_G))
error2 = np.sqrt(np.diag(pcov2))
params2 = unp.uarray(params2, error2)

x2 = np.linspace(0.36, 0.65, 1000)
# x2 = np.linspace(-0.62, -0.36, 1000)

m_gruen = unp.nominal_values(params2[0])
b_gruen = unp.nominal_values(params2[1])

fig3, ax3 = plt.subplots(layout='constrained')
ax3.plot(U_G, np.sqrt(I_G), 'x', color='green', label='Messdaten')
ax3.plot(x2, 
    linear(x2,m_gruen, b_gruen,),
    linestyle = '-',
    color = 'darkgreen',
    label = 'Lineare Regression'
)
ax3.legend(loc='best')
ax3.set_xlabel(r'$U\,$[V]')
ax3.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax3.grid()
ax3.set_ylim([-0.1e-6,4.1e-6])
ax3.set_xlim([0.36,0.63])
# ax3.set_xticks(np.arange(-0.62, -0.32, step=0.04))
fig3.savefig('Plots/gruen.pdf')

U_G_G = -params2[1]/ params2[0]
print(f'm und b für lin Reg Grün: {params2}')
print(f'Grenspannung Grün: U_G= {U_G_G}V')


# Orangenes Licht
U_O, I_O = np.genfromtxt('content/Messwerte/Orange.txt', unpack=True)
U_O = -1 * U_O

params3, pcov3 = curve_fit(linear, U_O, np.sqrt(I_O))
error3 = np.sqrt(np.diag(pcov3))
params3 = unp.uarray(params3, error3)

x3 = np.linspace(0.28, 0.55, 1000)
# x3 = np.linspace(-0.54, -0.28, 1000)

m_orange = unp.nominal_values(params3[0])
b_orange = unp.nominal_values(params3[1])

fig4, ax4 = plt.subplots(layout='constrained')
ax4.plot(U_O, np.sqrt(I_O), 'x', color='orangered', label='Messdaten')
ax4.plot(x3, 
    linear(x3, m_orange, b_orange,),
    linestyle = '-',
    color = 'orange',
    label = 'Lineare Regression'
)
ax4.legend(loc='best')
ax4.set_xlabel(r'$U\,$[V]')
ax4.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax4.grid()
ax4.set_ylim([-0.1e-6, 3.2e-6])
ax4.set_xlim([0.28, 0.55])
# ax4.set_xticks(np.arange(-0.54, -0.28, step=0.04))
fig4.savefig('Plots/orange.pdf')

U_G_O = -params3[1]/ params3[0]
print(f'm und b für lin Reg Orange: {params3}')
print(f'Grenspannung Orange: U_G= {U_G_O}V')


# Violettes Licht
U_V, I_V = np.genfromtxt('content/Messwerte/Violet.txt', unpack=True)
U_V = -1 * U_V

params4, pcov4 = curve_fit(linear, U_V, np.sqrt(I_V))
error4 = np.sqrt(np.diag(pcov4))
params4 = unp.uarray(params4, error4)

x4 = np.linspace(0.75, 1.22, 1000)
# x4 = np.linspace(-1.22, -0.75, 1000)

m_violet = unp.nominal_values(params4[0])
b_violet = unp.nominal_values(params4[1])

fig5, ax5 = plt.subplots(layout='constrained')
ax5.plot(U_V, np.sqrt(I_V), 'x', color='darkviolet', label='Messdaten')
ax5.plot(
    x4, 
    linear(x4, m_violet, b_violet,),
    linestyle = '-',
    color = 'indigo',
    label = 'Lineare Regression'
)
ax5.legend(loc='best')
ax5.set_xlabel(r'$U\,$[V]')
ax5.set_ylabel(r'$\sqrt{I}\,\left[\sqrt{\text{A}}\right]$')
ax5.grid()
ax5.set_ylim([-0.1e-6, 3.3e-6])
ax5.set_xlim([0.75, 1.22])
# ax5.set_xticks(np.arange(-1.22, -0.76, step=0.04))
fig5.savefig('Plots/violet.pdf')

U_G_V = -params4[1]/ params4[0]
print(f'm und b für lin Reg Violet: {params4}')
print(f'Grenspannung Violet: U_G= {U_G_V }V')



# Frequenzabhängigkeit
c = const.c
f_gruen = c/(546.07e-9)# nm
f_orange = c/576.96e-9 # nm
f_blau = c/435.83e-9 # nm
f_violet = c/407.78e-9 # nm 

f_all = np.array([f_gruen, f_orange, f_blau, f_violet])
U_G_all = np.array([U_G_G, U_G_O, U_G_B, U_G_V])
U_G_all = unp.nominal_values(U_G_all)

params5, pcov5 = curve_fit(linear, f_all, U_G_all)
error5 = np.sqrt(np.diag(pcov5))
params5 = unp.uarray(params5, error5)

x5 = np.linspace(0, 8e14, 1000)

m_planck = unp.nominal_values(params5[0])
b_planck = unp.nominal_values(params5[1])


fig6, ax6 = plt.subplots(layout = 'constrained')
ax6.plot(f_all, U_G_all, 'x', label = r'Berechnete $U_{\text{G}}$ der Lichtllinien')
ax6.plot(
    x5, 
    linear(x5, m_planck, b_planck),
    linestyle = '-',
    label = 'Lineare Regression'
)
ax6.legend(loc='best')
ax6.set_ylabel(r'$U_{\text{G}}\,$[eV]')
ax6.set_xlabel(r'$\nu\,$[Hz]')
ax6.grid()
ax6.set_xlim([0, 8e14])
ax6.set_ylim([-1.20, 1.5])
fig6.savefig('Plots/planck_berechnet.pdf')

print(f'm und b von f-U Diagramm berechnet {params5}')


f_all = np.array([f_gruen, f_orange, f_blau, f_violet])
U_G_all_gemessen = np.array([0.58, 0.50, 1.05, 1.1])

params6, pcov6 = curve_fit(linear, f_all, U_G_all_gemessen)
error6 = np.sqrt(np.diag(pcov6))
params6 = unp.uarray(params6, error6)


m_planck2 = unp.nominal_values(params6[0])
b_planck2 = unp.nominal_values(params6[1])

fig7, ax7 = plt.subplots(layout = 'constrained')
ax7.plot(f_all, U_G_all_gemessen, 'x', label = r'Gemessene $U_{\text{G}}$ der Lichtllinien')
ax7.plot(
    x5, 
    linear(x5, m_planck2, b_planck2),
    linestyle = '-',
    label = 'Lineare Regression'
)
ax7.legend(loc='best')
ax7.set_ylabel(r'$U_{\text{G}}\,$[eV]')
ax7.set_xlabel(r'$\nu\,$[Hz]')
ax7.grid()
ax7.set_xlim([0, 8e14])
ax7.set_ylim([-1.1, 1.4])
fig7.savefig('Plots/planck_gemessen.pdf')

print(f'm und b von f-U Diagramm gemessen {params6}')

h = const.h/const.e
phi_a = 4.05
phi_a2 = 4.6
phi_a_Ka = 2.25
print(f'e = {const.e}')
print(f'h = {h} eVs')
print(f'f g, o, b, v: {f_all}')


print('Diskussion --------------')

h_rel_gemessen = rel_Abweichung(params6[0], h)
h_rel_berechnet = rel_Abweichung(params5[0], h)
phi_rel_gemessen = rel_Abweichung(params6[1], phi_a)
phi_rel_berechnet = rel_Abweichung(params5[1], phi_a)
phi_rel_gemessen2 = rel_Abweichung(params6[1], phi_a2)
phi_rel_berechnet2 = rel_Abweichung(params5[1], phi_a2)
phi_rel_gemessen_Ka = rel_Abweichung(params6[1], phi_a_Ka)
phi_rel_berechnet_Ka = rel_Abweichung(params5[1], phi_a_Ka)

print(f'Abweichung h gemessen = {h_rel_gemessen}')
print(f'Abweichung h berechnet = {h_rel_berechnet}')
print(f'Abweichung Ag phi gemessen = {phi_rel_gemessen2} - {phi_rel_gemessen}')
print(f'Abweichung Ag phi berechnet = {phi_rel_berechnet2} - {phi_rel_berechnet}')
print(f'Abweichung Ka phi gemessen = {phi_rel_gemessen_Ka}')
print(f'Abweichung Ka phi berechnet = {phi_rel_berechnet_Ka}')
print(params6)
print(params5)