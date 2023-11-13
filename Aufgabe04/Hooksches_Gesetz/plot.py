import matplotlib.pyplot as plt
import numpy as np

x,F,D = np.genfromtxt("daten.txt", unpack=True)
params, covariance_matrix = np.polyfit(x, F, deg=1, cov=True)

print (params[0])
print (params[1])

x_plot = np.linspace(0,50)

fig, ax = plt.subplots(1,1, layout= "constrained")
ax.plot(x,F,"x", label="Messwerte")
ax.plot (
    x_plot,
    params[0] * x_plot + params[1],
    label = "Lineare Regression"
)
ax.set(
    xlabel = "Auslenkung der Feder $\Delta x \; [cm]$",
    ylabel = "Benötigte Kraft $F \; [N]$"
)
ax.legend();


fig.savefig("plot.pdf")