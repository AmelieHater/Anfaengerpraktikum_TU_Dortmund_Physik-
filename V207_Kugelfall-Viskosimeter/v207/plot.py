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
# r_kl = (15.57/2) * 10**(-1) # in cm (gemessen) +- 0,01mm
# r_gr = (15.76/2) * 10**(-1)# in cm (gemessen) +- 0.01mm
# print (f"Radius der kl. Kugel: {round(r_kl, 4)} und Durchmesser der kl. Kugel: {round(r_kl *2, 3)} (in cm)")
# print (f"Radius der gr. Kugel: {r_gr} und Durchmesser der gr. Kugel: {r_gr *2} (in cm)")

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

# Plot ln
temp, t_h1, t_r1, t_h2, t_r2, dichte = np.genfromtxt("temp.txt", unpack=True)
temp = unp.uarray(temp+273.15,1)
t_h = (t_h1 + t_h2) / 2
t_r = (t_r1 + t_r2) / 2
dichte_gr = np.ones(10) * dichte_gr
t_h_std = array (t_h1, t_h2)
print (t_h_std)
eta_hT = K_gr_h*(dichte_gr-dichte)*t_h
eta_rT = K_gr_r*(dichte_gr-dichte)*t_r



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