import matplotlib.pyplot as plt
import numpy as np

x,y = np.genfromtxt("content/Messwerte/Bandpass.txt", unpack = True)

fig, (ax1) = plt.subplots(1, 1, layout="constrained")
ax1.plot(x, y, "x",label="Kurve")
ax1.set_xlabel(1)
ax1.set_ylabel(2)
ax1.legend(loc="best")

fig.savefig("plot.pdf")