import homology2D as h2D
import homology3D as h3D
import numpy as np
import random as rd
from math import *

N =400


sphere = []
theta = 2 * np.pi * np.random.rand(N)
phi = np.arccos(1 - 2 * np.random.rand(N))

x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

for k in range(len(x)):
    p = h3D.Point(x[k], y[k], z[k], k)
    sphere.append(p)

torus = []

a = 1
b = 0.3

x = []
y = []
z = []

nbr_pts = 0

while nbr_pts < N:
    U = np.random.random()
    V = np.random.random()
    W = np.random.random()
    theta = 2 * np.pi * U
    phi = 2 * np.pi * V
    if W <=  (a+b*np.cos(theta))/(a+b):
        x.append((a + b * np.cos(theta)) * np.cos(phi))
        y.append((a + b * np.cos(theta)) * np.sin(phi))
        z.append(b * np.sin(theta))
        nbr_pts+=1

for k in range(len(x)):
    p = h3D.Point(x[k], y[k], z[k], k)
    torus.append(p)

# h3D.affiche(sphere, 'sphere.pdf')
#
# h3D.hom_pers(sphere, 0.10,'sphere_diag.pdf')

h3D.affiche(torus, 'torus.pdf')

h3D.hom_pers(torus, 0.05,'torus_diag.pdf',2)


# disque = []
# k = 0
# while k<N:
#     p = h2D.Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), k)
#     d = sqrt((p.x) ** 2 + (p.y) ** 2)
#     if (d < 0.8):
#         disque.append(p)
#         k += 1
# h2D.affiche(disque, 'disk.pdf')
# h2D.hom_pers(disque, 0.1,'disk_diag.pdf')
#
#
# anneau = h2D.annulus(N, 0.4, 0.6)
# h2D.affiche(anneau, 'anneau.pdf')
# h2D.hom_pers(anneau, 0.1,'anneau_diag.pdf')
#
#
# doubleanneau = []
# k = 0
# while k<N:
#     p = h2D.Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), k)
#     d1 = sqrt((p.x + 0.5) ** 2 + (p.y + 0.5) ** 2)
#     d2 = sqrt((p.x - 0.5) ** 2 + (p.y - 0.5) ** 2)
#     if (d1 > 0.4 and d1 < 0.5) or (d2 > 0.4 and d2 < 0.5):
#         doubleanneau.append(p)
#         k += 1
#
# h2D.affiche(doubleanneau, 'double_anneau.pdf')
# h2D.hom_pers(doubleanneau, 0.1,'double_anneau_diag.pdf')
