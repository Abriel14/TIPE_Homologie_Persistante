import Rips_complex as rc
import numpy as np
import numpy.random as rd

N = 20


theta = 2 * np.pi * rd.random(N)

r = np.sqrt(rd.random(N))

disk = []
for k in range(N):
    x = r[k] * np.cos(theta[k])
    y = r[k] * np.sin(theta[k])
    disk.append((x,y))

annulus = []
for k in range(N):
    x = (r[k]*0.2 + 1) * np.cos(theta[k])
    y = (r[k]*0.2 + 1) * np.sin(theta[k])
    annulus.append((x,y))

annulus_complex = rc.Rips_complex(annulus)

annulus_complex.show_pers_diagram()

# disk_complex = rc.Rips_complex(disk)
# disk_complex.show_pers_diagram()
