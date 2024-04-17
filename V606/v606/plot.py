import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

x, y = np.genfromtxt("content/Messwerte/Bandpass_weniger.txt", unpack=True)


# def gauss(x, a, v):
#    return np.exp(-a * (x-v)**2)
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(((x - mean) / 4 / stddev) ** 2))

#def gaussian(x, mean, stddev):
 #   return (1/stddev /np.sqrt(2 * np.pi)) * np.exp(-(((x - mean) / 4 / stddev) ** 2))


x_2 = np.linspace(26, 35, 100)
params, covariance_matrix = curve_fit(gaussian, x, y, p0=(130, 30,0.6))

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(x, y, "x", label="Kurve")
ax1.plot(x_2, gaussian(x_2, *params))
ax1.set_xlabel("f in kHz")
ax1.set_ylabel("U in Volt")
ax1.legend(loc="best")

fig.savefig("plot.pdf")
print(params)

# uncertainties = np.sqrt(np.diag(covariance_matrix))

# for name, value, uncertainty in zip("abc", params, uncertainties):
#    print(f"{name} = {value:8.3f} ± {uncertainty:.3f}")

# print(params)

# from scipy import optimize

# def gaussian(x, amplitude, mean, stddev):
#    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

# popt, _ = optimize.curve_fit(gaussian, x, y, p0=(0.5, 20))