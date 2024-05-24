import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

def rel_Abweichung(exp, theo):
    return (np.abs(exp-theo)/(theo)*100) # ist schon in Prozent

def fresnel_senkrecht(n, alpha):
    return (np.sqrt(n**2 - np.sin(alpha)**2)- np.cos(alpha))**2/(n**2 - 1)

def fresnel_parallel(n, alpha):
    return (n**2 * np.cos(alpha) - np.sqrt(n**2 - np.sin(alpha)**2))/(n**2 * np.cos(alpha) + np.sqrt(n**2 - np.sin(alpha)**2))

def n_senkrecht(E,alpha):
    return np.sqrt(np.abs((E**2 + 2*E*np.cos(2*alpha) + 1)/(E**2 - 2*E + 1)))

def n_parallel(E, alpha):
    return (((E+1)/(E-1))**2 * (1/(2* np.cos(alpha)**2)) + np.sqrt(np.abs((1/(4*np.cos(alpha)**2)) * ((E+1)/(E-1))**4 - ((E+1)/(E-1))**2 * np.tan(alpha)**2)))

winkel, I_senkrecht, I_parallel = np.genfromtxt('content/Messwerte/Messwerte.txt', unpack = True)
winkel = winkel * np.pi/180
I_dunkel = 2.8 * 10**(-9) # Ampere
I_0 = 0.46 * 10**(-3) # Ampere
alpha_B = 75 * np.pi/180

n_theo = np.tan(alpha_B)
print(f'I_0: {I_0/10**-6}')
print(f'Brewster: {alpha_B}')
print(f'n_B_theo: {n_theo}')
I_plot_s = np.sqrt(I_senkrecht / I_0)
I_plot_p = np.sqrt(I_parallel / I_0)
# print(f'I_p: {I_plot_p}')

n_s = n_senkrecht(I_plot_s, winkel)
n_p = n_parallel(I_plot_p, winkel)
# print(f'n_s: {n_s}')
# print(f'n_p: {n_p}')

n_s_mean = ufloat(np.mean(n_s), np.std(n_s))
n_p_mean = ufloat(np.mean(n_p[n_p < 10 ]), np.std(n_p[n_p < 10]))
print(f'n_s mean: {n_s_mean}')
print(f'n_p mean: {n_p_mean}')

x = np.linspace(0, 90, 1000)
x = x * np.pi/180

x_1 = np.linspace(0, alpha_B, 1000)
x_2 = np.linspace(alpha_B, 87, 2000)
fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(winkel * 180/np.pi, I_plot_s, 'x', color='cornflowerblue', label='Messdaten senkrecht')
ax1.plot(x * 180/np.pi, fresnel_senkrecht(n_theo, x), '-', color='indigo', label='Theoriekurve senkrecht')
ax1.plot(winkel * 180/np.pi, I_plot_p, 'x', color='darkorange', label='Messdaten parallel')
ax1.plot(x_1 * 180/np.pi, fresnel_parallel(n_theo , x_1), '-', color='orangered', label='Theoriekurve parallel')
ax1.plot(x_2 * 180/np.pi, -1*fresnel_parallel(n_theo, x_2), '-', color='orangered')
ax1.vlines(75, ymin = 0, ymax = 1.05, linestyle = ':', label='Brewsterwinkel')
ax1.set_xlabel(r'Winkel $\alpha\,[°]$')
ax1.set_ylabel(r'Intensität $I\,[\text{A}]$')
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax1.grid()
ax1.legend(loc='best')
ax1.set_xlim([0,90])
ax1.set_ylim([-0.05,1.05])
fig.savefig('plot.pdf')

print(f'Diskussion------------')
n_rel_s = rel_Abweichung(np.mean(n_s), n_theo)
n_rel_p_reduziert = rel_Abweichung(np.mean(n_p[n_p < 10]), n_theo)
n_rel_p = rel_Abweichung(np.mean(n_p), n_theo)
print(f'rel Abweichugn n_s: {n_rel_s}')
print(f'rel Abweichugn n_p reduziert: {n_rel_p_reduziert}')
print(f'rel Abweichugn n_p: {n_rel_p}')