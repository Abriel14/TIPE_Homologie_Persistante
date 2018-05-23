import homology as h
import numpy.random as rd
import math
N = 600

h.limite(100,250,0.5,0.7,10)









# points = [h.Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), i) for i in range(N)]
#
# anneau = h.annulus(N,0.5,0.6)
#
# h.affiche(anneau)
#
#
#
# cplx,nbL = h.Lazy_Witness_Complex(anneau , 0.1, 1)
# D = h.calcule_D(cplx,nbL)
# DD = h.reduction_D(D)
# C = h.paires_pers(D,cplx,nbL)
# h.diag_pers(C)
#
#
# doubleanneau=[]
#
# k= 0
# for i in range(2*N) :
#     p = h.Point(-1 + 2*rd.random(),-1 + 2*rd.random(),k)
#     d1 = math.sqrt((p.x+0.5)**2+(p.y+0.5)**2)
#     d2 = math.sqrt((p.x-0.5)**2+(p.y-0.5)**2)
#     if (d1 > 0.4 and d1 < 0.5) or (d2 > 0.4 and d2 < 0.5) :
#         doubleanneau.append(p)
#         k+=1
#
# h.affiche(doubleanneau)
#
# cplx,nbL = h.Lazy_Witness_Complex(doubleanneau , 0.1, 1)
# D = h.calcule_D(cplx,nbL)
# DD = h.reduction_D(D)
# C = h.paires_pers(D,cplx,nbL)
# h.diag_pers(C)