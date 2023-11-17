import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp


# Messdaten
m_kl = 4.4531 # in g (gegeben)
m_gr = 4.9528 # in g (gegeben)

d_kl = ufloat(15.57,0.01) * 10**-1 # in cm (gemessen)
d_gr = ufloat(15.76,0.01) * 10**-1 # in cm (gemessen)
print (f"Durchmesser kl. Kugel: {d_kl} (in cm)")
print (f"Durchmesser gr. Kugel: {d_gr} (in cm)")

# Dichte \rho = m/V 
V_kl = (4/3) * np.pi * (d_kl/2)**3
V_gr = (4/3) * np.pi * (d_gr/2)**3
dichte_kl = m_kl/(V_kl)
dichte_gr = m_gr/(V_gr)
print (f"Volumen der kl. Kugel: {V_kl} (in cm^3)")
print (f"Volumen der gr. Kugel: {V_gr} (in cm^3)")
print (f"Dichte der kl. Kugel: {dichte_kl} (in g/cm^3)")
print (f"Dichte der gr. Kugel: {dichte_gr} (in g/cm^3)")

# Fallzeit kl. Kugel
t_kl_r0, t_kl_h0 = np.genfromtxt("kl_Kugel.txt", unpack=True)
t_gr_h0, t_gr_r0 = np.genfromtxt ("gr_Kugel.txt", unpack=True)
t_kl_r_mean = np.mean(t_kl_r0)
t_kl_h_mean = np.mean(t_kl_h0)
t_gr_r_mean = np.mean(t_gr_r0)
t_gr_h_mean = np.mean(t_gr_h0)
t_kl_r_std = np.std(t_kl_r0)
t_kl_h_std = np.std(t_kl_h0)
t_gr_r_std = np.std(t_gr_r0)
t_gr_h_std = np.std(t_gr_h0)

t_kl_r = ufloat(t_kl_r_mean, t_kl_r_std)
t_kl_h = ufloat(t_kl_h_mean, t_kl_h_std)
t_gr_r = ufloat(t_gr_r_mean, t_gr_r_std)
t_gr_h = ufloat(t_gr_h_mean, t_gr_h_std)

print (f"Gemittelte Fallzeit (hoch) kl: {t_kl_h}, gr: {t_gr_h}")
print (f"Gemittelte Fallzeit (runter) kl: {t_kl_r}, gr: {t_gr_r}")

# Apparaturkonstanten
K_kl = 0.07640 # in m*Pa*cm^3/g (gegeben)
dichte_wasser = 0.998207 # in g/cm^3 (Internet)
eta_h = K_kl * (dichte_kl-dichte_wasser) * t_kl_h
eta_r =K_kl * (dichte_kl-dichte_wasser) * t_kl_r

K_gr_h = eta_h/((dichte_gr-dichte_wasser) * t_gr_h)
K_gr_r = eta_r/((dichte_gr-dichte_wasser) * t_gr_r)

print(f"Viskosität hoch: {eta_h} (in mPa*s)")
print (f"Viskosität runter: {eta_r} (in mPa*s)")
print (f"Apparaturkonstante K_gr_h: {K_gr_h} (in mPa*cm^3/g)")
print (f"Apparaturkonstante K_gr_r: {K_gr_r} (in mPa*cm^3/g)")

# Reynoldsche Zahl
Re_kl_h = 100*(dichte_wasser * (10/ t_kl_h) * d_gr) / eta_h
Re_kl_r = 100*(dichte_wasser * (10/ t_kl_r) * d_gr) / eta_r
Re_gr_h = 100*(dichte_wasser * (5/ t_gr_h) * d_gr) / eta_h
Re_gr_r = 100*(dichte_wasser * (5/ t_gr_r) * d_gr) / eta_r

print (f"Reynoldsche Zahl Re_kl_h: {Re_kl_h}")
print (f"Reynoldsche Zahl Re_kl_r: {Re_kl_r}")
print (f"Reynoldsche Zahl Re_gr_h: {Re_gr_h}")
print (f"Reynoldsche Zahl Re_gr_r: {Re_gr_r}")

# Temperaturabhängigkeit
temp, t_h1, t_r1, t_h2, t_r2, dichte = np.genfromtxt("temp.txt", unpack=True)
temp = unp.uarray(temp+273.15,1)
dichte_gr = np.ones((1,10))* dichte_gr

t_ges = np.matrix((t_h1 , t_r1, t_h1, t_r1))
print (f"t_ges: {t_ges}")
t_ges_st = np.std (t_ges, axis=0)
t_ges_mean = np.mean (t_ges, axis=0)

t_ges = unp.uarray(t_ges_mean, t_ges_st)

# Viskosität abhängig Temp.
dichten = (dichte_gr -dichte).reshape (10,1)

eta_T = K_gr_h * (dichten) * t_ges

eta_T = np.diagonal(eta_T)

print (f"eta_T {eta_T}")

# Plots
x_err = 1/temp
y_err = unp.log(eta_T)
x = 1/unp.nominal_values(temp)
y = np.log(unp.nominal_values(eta_T))

print (f"temp: {temp}")
print (f"y: {y}")
params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
#params, covariance_matrix = np.polyfit(x, y, deg=1, cov=True)
errors = np.sqrt(np.diag(covariance_matrix))

for name, value, error in zip ('ab', params, errors):
    print(f"{name} = {value:.3f} ± {error:.3f}")

x_plot = np.linspace(31e-4,33e-4,10000)
# fig, ax = plt.subplots(1,1, layout="constrained")

#plt.plot (x, y, "x", label = "Messwerte")
plt.errorbar(x,y, xerr=unp.std_devs(x_err), yerr=unp.std_devs(y_err), fmt="x", label = "Messwerte")
plt.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label ="Lineare Regression",
    linewidth=1,
)
plt.grid()
plt.xlabel(r'$T^{-1}$ [K $^{-1}$]')
plt.ylabel(r"$\ln{ \left ( \eta \right )}$ [mPa$\cdot$s]")
plt.legend(loc = "best")
plt.margins(0.075)
plt.savefig("plot.pdf")

a = ufloat(1680.367, 30.445)
b = ufloat (-5.567, 0.097)

print (f"B = {a}")
print (f"A = {unp.exp(b)}")
# t_ges_st = np.std (t_ges, )
# t_h0 = np.matrix((t_h1, t_h2))
# t_h_st = np.std(t_h0, axis=0)
# t_h_m = np.mean(t_h0, axis = 0)

# t_r0 = np.matrix ((t_r1, t_r2))
# t_r_st = np.std(t_r0, axis = 0)
# t_r_m = np.mean(t_r0, axis = 0)
# # print (f"std hoch : {t_h_st}")
# # print (f"std runter: {t_r_st}")
# t_h = unp.uarray(t_h_m, t_h_st)
# t_r = unp.uarray(t_r_m, t_r_st)
# # print (f"mean hoch : {t_h_m}")
# # print (f"mean runter: {t_r_m}")

# print (f"t hoch : {t_h}")
# print (f"t runter: {t_r}")
# dichte_gr = np.ones((1,10))* dichte_gr
# print (f"dichte: {dichte_gr}")
# eta_hT = np.zeros(10)
# eta_rT = np.zeros(10)
# print(f"dichte_gr: {dichte_gr}")
# print(f"dichte: {dichte}")
# print (f"t_h: {t_h}")
# for x in range(10):
#    print(f"x: {x}")
#    print(dichte_gr)
#    print(K_gr_h)
#    print(K_gr_r)
#    print(eta_hT[x])
#    print(dichte[x])
#    print (f"t_h:{t_h[0,x]}")
#    eta_hT[x] = K_gr_h * (dichte_gr - dichte[x]) * t_h[0,x]
#    eta_rT[x] = K_gr_r * (dichte_gr - dichte[x]) * t_r[0,x]

# x = dichte_gr - dichte
# print (f"x: {x}")
# y= t_h * (x)
# print (y)
# eta_hT = K_gr_h * t_h * (x)
# eta_rT = K_gr_r * t_r * (x)



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