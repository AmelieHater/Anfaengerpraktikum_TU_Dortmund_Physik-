import numpy as np
import scipy.constants
from uncertainties import ufloat
import uncertainties.unumpy as unp


sigma = scipy.constants.sigma
h = scipy.constants.h
k = scipy.constants.k
m0 = scipy.constants.m_e
pi = scipy.constants.pi
e0 = scipy.constants.elementary_charge
eta = 0.28
f = 3.2e-5
#b = ufloat(7.0364, 0.3406)
N_WL = 0.95
print("sigma: ", sigma)
print("h: ", h)
print("k: ", k)
print("m0: ", m0)
print("pi: ", pi)
print("e0: ", e0)
print("eta: ", eta)
print("f: ", f)

I_s = np.array([60e-6, 140e-6, 319e-6, 719e-6, 1356e-6])

I_h_1 = np.array([2, 2.1, 2.2, 2.3, 2.4])
I_h_err = np.full_like(I_h_1, 0.1)
I_h = unp.uarray(I_h_1, I_h_err)

U_h_1 = np.array([4, 4, 4.5, 5, 5])
U_h_err = np.full_like(U_h_1, 0.5)
U_h = unp.uarray(U_h_1, U_h_err)


T_L = ((U_h*I_h-N_WL)/(f*eta*sigma))**0.25
W = -k*T_L*(unp.log(I_s*h**3/(4*pi*e0*m0*k**2*f*T_L**2)))/e0
W_mean = W.mean()
#T_k = e0/(k*b)


print(f'Kathodentemperatur über Leistungsbilanz in Kelvin: \n {T_L.reshape(-1,1)}')
#print(f'Kathodentemperatur über Anlaufstromgebiet: {T_k} K')
print(f'Austrittsarbeiten in eV: \n {W}') 
print(f'Mittelwert: {W_mean} eV')




