import Rips_complex as rc
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
N = 10



theta = 2 * np.pi * rd.random(N)

r = np.sqrt(rd.random(N))

disk = []
x_plot = []
y_plot = []


for k in range(N):
    x = r[k] * np.cos(theta[k])
    y = r[k] * np.sin(theta[k])
    x_plot.append(x)
    y_plot.append(y)
    disk.append((x,y))


annulus = []

for k in range(N):
    x = (r[k]*0.01 + 0.5) * np.cos(theta[k])
    y = (r[k]*0.01 + 0.5) * np.sin(theta[k])
    annulus.append((x,y))


plt.plot(x_plot,y_plot,'o')
plt.show()

annulus_complex = rc.Rips_complex(annulus)
annulus_complex.show_pers_diagram()

disk_complex = rc.Rips_complex(disk)
disk_complex.show_pers_diagram()
