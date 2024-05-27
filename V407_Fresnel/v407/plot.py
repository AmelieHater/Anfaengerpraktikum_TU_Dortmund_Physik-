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
    return np.sqrt(((E**2 - 2*E*np.cos(2*alpha) + 1)/(E**2- 2*E + 1)))

# def n_parallel(E, alpha):
#     return np.sqrt(((1+E)/(1-E))**2 * (1/(2* np.cos(alpha)**2)) + (1/(2*np.cos(alpha)**2 *(1-E)**2))*np.sqrt((1+E)**2-4*np.cos(alpha)**2*(1-E)**2*np.sin(alpha)**2))

def n_parallel(E,alpha):
    return np.sqrt((1+E)**2/(1-E)**2 * 1/(2*np.cos(alpha)**2) + np.sqrt((1+E)**2/(4*np.cos(alpha)**4*(1-E)**4)-1/(1-E)**2 * np.tan(alpha)**2))

winkel, I_senkrecht, I_parallel = np.genfromtxt('content/Messwerte/Messwerte.txt', unpack = True)
winkel = winkel * np.pi/180
I_dunkel = 2.8 * 10**(-9) # Ampere
I_0 = 0.46 * 10**(-3) # Ampere
alpha_B = 75 * np.pi/180

n_theo = np.tan(alpha_B)
print(f'I_0: {I_0/10**-6}')
print(f'Brewster: {alpha_B * 180/np.pi}')
print(f'n_B_theo: {n_theo}')
I_plot_s = np.sqrt(I_senkrecht / I_0)
I_plot_p = np.sqrt(I_parallel / I_0)
# print(f'I_p: {I_plot_p}')

n_s = n_senkrecht(I_plot_s, winkel)
n_p = n_parallel(I_plot_p, winkel)
print(f'n_s: {np.round(n_s,2)}')
print(f'n_p: {np.round(n_p,2)}')

n_s_mean = ufloat(np.mean(n_s[n_s < 4]), np.std(n_s[n_s < 4]))
n_p_mean = ufloat(np.mean(n_p[n_p < 4]), np.std(n_p[n_p < 4]))
n_s_mean_schlecht = ufloat(np.mean(n_s), np.std(n_s))
n_p_mean_schlecht = ufloat(np.mean(n_p), np.std(n_p))

print(f'n_s mean: {n_s_mean}')
print(f'n_p mean: {n_p_mean}')
print(f'n_p mean schlecht: {n_p_mean_schlecht}')
print(f'n_s mean schlecht: {n_s_mean_schlecht}')

x = np.linspace(0, 90, 1000)
x = x * np.pi/180

x_1 = np.linspace(0, 67.7, 1000)
x_2 = np.linspace(67.7, 90, 2000)
fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(winkel * 180/np.pi, I_plot_s, 'x', color='cornflowerblue', label='Messdaten senkrecht')
ax1.plot(x * 180/np.pi, fresnel_senkrecht(unp.nominal_values(n_s_mean), x), '-', color='indigo', label='Theoriekurve senkrecht')
ax1.plot(winkel * 180/np.pi, I_plot_p, 'x', color='darkorange', label='Messdaten parallel')
# ax1.plot(x * 180/np.pi, fresnel_parallel(np.mean(n_p[n_p < 6 ]) , x), '-', color='orangered', label='Theoriekurve parallel')
ax1.plot(x_1, fresnel_parallel(unp.nominal_values(n_p_mean) , x_1*np.pi/180), '-', color='orangered', label='Theoriekurve parallel')
ax1.plot(x_2, -1*fresnel_parallel(unp.nominal_values(n_p_mean), x_2*np.pi/180), '-', color='orangered')
ax1.vlines(75, ymin = 0, ymax = 1.05, linestyle = ':', label='Brewsterwinkel')
ax1.set_xlabel(r'$\alpha\,[°]$')
ax1.set_ylabel(r'$\sqrt{\frac{I_{r}(\alpha)}{I_0}}$')
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax1.grid()
ax1.legend(loc='best')
ax1.set_xlim([0,90])
ax1.set_ylim([-0.05,1.05])
fig.savefig('plot.pdf')

print(f'Diskussion------------')
n_rel_s_reduziert = rel_Abweichung(unp.nominal_values(n_s_mean), n_theo)
n_rel_p_reduziert = rel_Abweichung(unp.nominal_values(n_p_mean), n_theo)
n_rel_p = rel_Abweichung(unp.nominal_values(n_p_mean_schlecht), n_theo)
n_rel_s = rel_Abweichung(unp.nominal_values(n_s_mean_schlecht), n_theo)
print(f'rel Abweichugn n_s reduziert: {n_rel_s_reduziert}')
print(f'rel Abweichugn n_p reduziert: {n_rel_p_reduziert}')
print(f'rel Abweichugn n_p: {n_rel_p}')
print(f'rel Abweichugn n_s: {n_rel_s}')