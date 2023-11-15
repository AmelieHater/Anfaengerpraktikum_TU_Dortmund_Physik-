import matplotlib.pyplot as plt
import numpy as np

# Messdaten
m_kl = 4.4531 # in g (gegeben)
m_gr = 4.9528 # in g (gegeben)

r_kl = (15.57/2) * 10**(-1) # in cm (gemessen) +- 0,01mm
r_gr = (15.76/2) * 10**(-1)# in cm (gemessen) +- 0.01mm
print (f"Radius der kl. Kugel: {round(r_kl, 4)} und Durchmesser der kl. Kugel: {round(r_kl *2, 3)} (in cm)")
print (f"Radius der gr. Kugel: {r_gr} und Durchmesser der gr. Kugel: {r_gr *2} (in cm)")

# Dichte \rho = m/V mit V=4/3*pi*r^3
dichte_kl = round(m_kl/((4/3) * np.pi * (r_kl)**3), 6)
dichte_gr = round(m_gr/((4/3) * np.pi * (r_gr)**3), 6)
print (f"Dichte der kl. Kugel: {dichte_kl} (in g/cm^3)")
print (f"Dichte der gr. Kugel: {dichte_gr} (in g/cm^3)")

# Fallzeit kl. Kugel
t_kl_r, t_kl_h = np.genfromtxt("kl_Kugel.txt", unpack=True)
t_gr_h, t_gr_r = np.genfromtxt ("gr_Kugel.txt", unpack=True)
t_kl_r_mittel = sum(t_kl_r)/len(t_kl_r)
t_kl_h_mittel = sum(t_kl_h)/len(t_kl_h)
t_gr_r_mittel = sum(t_gr_r)/len(t_gr_r)
t_gr_h_mittel = sum(t_gr_h)/len(t_gr_h)

print (f"Gemittelte Fallzeit (hoch) kl: {t_kl_h_mittel}, gr: {t_gr_h_mittel}")
print (f"Gemittelte Fallzeit (runter) kl: {t_kl_r_mittel}, gr: {t_gr_r_mittel}")

# Apparatekonstanten
K_kl = 0.07640 # in m*Pa*cm^3/g (gegeben)
dichte_wasser = 0.998207 # in g/cm^3 (Internet)
eta_h = K_kl * (dichte_kl-dichte_wasser) * t_kl_h_mittel
eta_r =K_kl * (dichte_kl-dichte_wasser) * t_kl_r_mittel

K_gr_h = eta_h/((dichte_gr-dichte_wasser) * t_gr_h_mittel)
K_gr_r = eta_r/((dichte_gr-dichte_wasser) * t_gr_r_mittel)

print(f"eta_hoch: {eta_h} (in mPa*s)")
print (f"eta_runter: {eta_r} (in mPa*s)")
print (f"K_gr_h: {K_gr_h} (in mPa*cm^3/g)")
print (f"K_gr_r: {K_gr_r} (in mPa*cm^3/g)")

# Reynoldsche Zahl
Re_kl_h = 100*(dichte_wasser * (10/t_kl_h_mittel) * r_gr) / eta_h
Re_kl_r = 100*(dichte_wasser * (10/ t_kl_r_mittel) * r_gr) / eta_r
Re_gr_h = 100*(dichte_wasser * (5/ t_gr_h_mittel) * r_gr) / eta_h
Re_gr_r = 100*(dichte_wasser * (5/ t_gr_r_mittel) * r_gr) / eta_r

print (f"Re_kl_h: {Re_kl_h}")
print (f"Re_kl_r: {Re_kl_r}")
print (f"Re_gr_h: {Re_gr_h}")
print (f"Re_gr_r: {Re_gr_r}")

# Plot ln

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