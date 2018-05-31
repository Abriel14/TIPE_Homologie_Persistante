import homology2D as h2D
import homology3D as h3D
import numpy as np
import math
import random as rd

# from test import reduction_D

N = 600

# h.limite(100,250,0.5,0.7,10)


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

a = 0.7
b = 0.3

theta = 2 * np.pi * np.random.rand(N)
phi = 2 * np.pi * np.random.rand(N)

x = (a + b * np.cos(phi)) * np.cos(theta)
y = (a + b * np.cos(phi)) * np.sin(theta)
z = b * np.sin(phi)

for k in range(len(x)):
    p = h3D.Point(x[k], y[k], z[k], k)
    torus.append(p)

# h3D.affiche(sphere, 'sphere.pdf')
#
# cplx, nbL = h3D.Lazy_Witness_Complex(sphere, 0.05, 1)
# D = h3D.calcule_D(cplx, nbL)
# DD = h3D.reduction_D(D)
# C = h3D.paires_pers(D, cplx, nbL)
# h3D.diag_pers(C, 'sphere_diag.pdf')
#
# h3D.affiche(torus, 'torus.pdf')
#
# cplx, nbL = h3D.Lazy_Witness_Complex(torus, 0.05, 4)
# D = h3D.calcule_D(cplx, nbL)
# DD = h3D.reduction_D(D)
# C = h3D.paires_pers(D, cplx, nbL)
# h3D.diag_pers(C, 'torus_diag.pdf')



disque = []
k = 0
for i in range(N):
    p = h2D.Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), k)
    d = math.sqrt((p.x) ** 2 + (p.y) ** 2)
    if (d < 0.8):
        disque.append(p)
        k += 1
h2D.affiche(disque, 'disk.pdf')

cplx, nbL = h2D.Lazy_Witness_Complex(disque, 0.1, 1)
D = h2D.calcule_D(cplx, nbL)
DD = h2D.reduction_D(D)
C = h2D.paires_pers(D, cplx, nbL)
h2D.diag_pers(C, 'disk_diag.pdf')

anneau = h2D.annulus(N, 0.5, 0.6)

h2D.affiche(anneau, 'anneau.pdf')

cplx, nbL = h2D.Lazy_Witness_Complex(anneau, 0.1, 1)
D = h2D.calcule_D(cplx, nbL)
DD = h2D.reduction_D(D)
C = h2D.paires_pers(D, cplx, nbL)
h2D.diag_pers(C, 'anneau_diag.pdf')

doubleanneau = []

k = 0
for i in range(2 * N):
    p = h2D.Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), k)
    d1 = math.sqrt((p.x + 0.5) ** 2 + (p.y + 0.5) ** 2)
    d2 = math.sqrt((p.x - 0.5) ** 2 + (p.y - 0.5) ** 2)
    if (d1 > 0.4 and d1 < 0.5) or (d2 > 0.4 and d2 < 0.5):
        doubleanneau.append(p)
        k += 1

h2D.affiche(doubleanneau, 'double_anneau.pdf')

cplx, nbL = h2D.Lazy_Witness_Complex(doubleanneau, 0.1, 1)
D = h2D.calcule_D(cplx, nbL)
DD = h2D.reduction_D(D)
C = h2D.paires_pers(D, cplx, nbL)
h2D.diag_pers(C, 'double_anneau_diag.pdf')
