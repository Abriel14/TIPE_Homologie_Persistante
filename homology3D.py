import numpy as np
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import math as m
from operator import itemgetter
import random as rd
import math


class Point:
    """ Point class represents and manipulates x,y coords. """

    def __init__(self, x, y, z, numero):
        """ Create a new point en x,y """
        self.x = x
        self.y = y
        self.z = z
        self.numero = numero

    def dist(self, p):
        """calcule la distance entre un point et un autre"""
        d = m.sqrt((self.x - p.x) ** 2 + (self.y - p.y) ** 2 + (self.z - p.z) ** 2)
        return d


def affiche(tabP, scatt=20):
    """affiche un nuage de point"""
    px = [tabP[i].x for i in range(len(tabP))]
    py = [tabP[i].y for i in range(len(tabP))]
    pz = [tabP[i].z for i in range(len(tabP))]
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.scatter(px, py, pz)
    show()


def mini_maxi(tabP, delta, i, j, nbrL):
    """calcule le minimum des max(d(a1,w),d(a2,w)) - delta[w] pour tous les point w n'appartenants pas à la filtration"""
    res = [max(tabP[w].dist(tabP[i]), tabP[w].dist(tabP[j])) - delta[w - nbrL] for w in range(nbrL, len(tabP))]
    return (min(res))


def Lazy_Witness_Complex(tabP, pourcentage, v):
    """crée le Lazy witness Complex"""
    nbrL = int(len(tabP) * pourcentage)
    # représente le nombre de points du filtrage
    D = [[tabP[k].dist(tabP[j]) for k in range(nbrL)] for j in range(nbrL, len(tabP))]
    # D[i,j] représente d(ai,aj) pour ai dans le filtrage et ai n'étant pas dans le filtrage
    for i in range(len(D)):
        D[i].sort()
    delta = [D[k][v] for k in range(len(D))]
    # delta[w,v] représente la v-ème plus petite distance entre le point témoin w n'étant pas dans le filtrage et le filtrage
    date_apparition = [[mini_maxi(tabP, delta, i, j, nbrL) for i in range(nbrL)] for j in range(nbrL)]
    # Tableau de la date d'apparition des 1_simplexes
    simplexes = []
    for i in range(nbrL):
        for j in range(i + 1, nbrL):
            simplexes.append([(i, j), date_apparition[i][j], 1])
    for i in range(nbrL):
        for j in range(i + 1, nbrL):
            for k in range(j + 1, nbrL):
                simplexes.append(
                    [(i, j, k), max(date_apparition[i][j], date_apparition[i][k], date_apparition[j][k]), 2])
    for i in range(nbrL):
        for j in range(i + 1, nbrL):
            for k in range(j + 1, nbrL):
                for l in range(k + 1, nbrL):
                    simplexes.append(
                        [(i, j, k, l),
                         max(date_apparition[i][j], date_apparition[i][k], date_apparition[i][l], date_apparition[j][k],
                             date_apparition[j][l], date_apparition[k][l]), 3])

    # simplexes est un grand tableau regroupant les informations [simplexe, date, type de simplexe] sur les simplexes créés
    simplexes.sort(key=itemgetter(1))
    # on le trie par date d'apparition des simplexes croissante
    return (simplexes, nbrL)


def calcule_D(lazycplx, nbrL):
    """calcule la matrice D"""
    D = np.zeros((len(lazycplx) + nbrL, len(lazycplx) + nbrL))
    pos_unsimplexe = np.zeros((nbrL, nbrL))
    pos_deuxsimplexe = np.zeros((nbrL, nbrL, nbrL))
    for v in range(len(lazycplx)):
        pos = v + nbrL
        (splx, time, type) = lazycplx[v]
        if type == 1:
            # c'est un 1 simplexe
            pos_unsimplexe[splx[0], splx[1]] = pos
            pos_unsimplexe[splx[1], splx[0]] = pos
            D[splx[0], pos] = 1
            D[splx[1], pos] = 1
        elif type == 2:
            # c'est un 2 simplexe
            pos_deuxsimplexe[splx[0], splx[1], splx[2]] = pos
            pos_deuxsimplexe[splx[0], splx[2], splx[1]] = pos
            pos_deuxsimplexe[splx[1], splx[0], splx[2]] = pos
            pos_deuxsimplexe[splx[1], splx[2], splx[0]] = pos
            pos_deuxsimplexe[splx[2], splx[1], splx[0]] = pos
            pos_deuxsimplexe[splx[2], splx[0], splx[1]] = pos
            D[int(pos_unsimplexe[splx[0], splx[1]]), pos] = 1
            D[int(pos_unsimplexe[splx[1], splx[2]]), pos] = 1
            D[int(pos_unsimplexe[splx[0], splx[2]]), pos] = 1
        else:
            # c'est un 3 simplexe
            D[int(pos_deuxsimplexe[splx[0], splx[1], splx[2]]), pos] = 1
            D[int(pos_deuxsimplexe[splx[0], splx[1], splx[3]]), pos] = 1
            D[int(pos_deuxsimplexe[splx[0], splx[2], splx[3]]), pos] = 1
            D[int(pos_deuxsimplexe[splx[1], splx[2], splx[3]]), pos] = 1
    return (D)


def low(mat, j):
    """calcule low(j), low(j) est le plus grand indice i tel que la i-ème valaur de la colonne j de mat vaut 1 """
    res = -1
    for i in range(len(mat)):
        if mat[i][j] == 1:
            res = i
    return (res)


def reduction_D(mat):
    """procède à la réduction de la matrice D"""
    listlow = [low(mat, j) for j in range(len(mat))]
    for j in range(len(mat)):
        ilexiste = True
        while listlow[j] != -1 and ilexiste:
            ilexiste = False
            j0 = 0
            while j0 < j and not (ilexiste):
                if listlow[j0] == listlow[j]:
                    ilexiste = True
                    mat[:, j] = (mat[:, j] + mat[:, j0]) % 2
                    listlow[j] = low(mat, j)
                else:
                    j0 = j0 + 1
    return (mat)


def paires_pers(D, cplx, nbrL):
    listlow = [low(D, j) for j in range(len(D))]
    res = []
    for i in range(len(listlow)):
        if listlow[i] != -1:
            (splx, time_d, type) = cplx[i - nbrL]
            (splx, time_b, typeb) = cplx[listlow[i] - nbrL]
            res.append(((time_b, time_d), type))
    return (res)


def diag_pers(paires):
    h0_birth = []
    h0_death = []
    h1_birth = []
    h1_death = []
    h2_birth = []
    h2_death = []
    for p in paires:
        [(b, d), type] = p
        if type == 1:
            h0_birth.append(b)
            h0_death.append(d)
        if type == 2:
            h1_birth.append(b)
            h1_death.append(d)
        if type == 3:
            h2_birth.append(b)
            h2_death.append(d)
    figure()
    subplot(311)
    title("H0")
    xlabel('Indice du simplexe')
    ylabel("Indice de la classe d'homologie")
    hlines(range(len(h0_birth)), h0_birth, h0_death, colors='b')
    subplot(312)
    title("H1")
    xlabel('Indice du simplexe')
    ylabel("Indice de la classe d'homologie")
    hlines(range(len(h1_birth)), h1_birth, h1_death, colors='b')
    subplot(313)
    title("H2")
    xlabel('Indice du simplexe')
    ylabel("Indice de la classe d'homologie")
    hlines(range(len(h2_birth)), h2_birth, h2_death, colors='b')
    show()



# def calcul_D(lzyCplx, nbr_tem):
#     D = [[] for k in range(nbr_tem + len(lzyCplx))]
#     un_splx_pos = np.zeros((nbr_tem, nbr_tem))
#     for k in range(len(lzyCplx)):
#         position = nbr_tem + k
#         data, time, type = lzyCplx[k]
#         if type == 1:
#             un_splx_pos[data[0], data[1]] = position
#             un_splx_pos[data[1], data[0]] = position
#             D[position] = sorted([data[0], data[1]])
#         if type == 2:
#             D[position] = sorted(
#                 [un_splx_pos[data[0], data[1]], un_splx_pos[data[0], data[2]], un_splx_pos[data[1], data[2]]])
#     return (D)
#
#
# def low_j(M, j):
#     if M[j] != []:
#         return (sorted(M[j])[-1])
#     else:
#         return (0)
#
#
# def sum_mod_2(M, i, j):
#     column_i = sorted(M[i])
#     column_j = sorted(M[j])
# #    print(column_i,column_j)
#     if len(column_j) < len(column_i):
#         temp = column_j
#         column_j = column_i
#         column_i = temp
#     sum = []
#     while column_j != []:
#         if column_i != []:
#             a = column_i[-1]
#             b = column_j[-1]
#             if a < b:
#                 sum.append(column_j.pop())
#             if a > b:
#                 sum.append(column_i.pop())
#             if a == b:
#                 column_i.pop()
#                 column_j.pop()
#         else:
#             sum.append(column_j.pop())
# #    print(sorted(sum))
#     return (sorted(sum))
#
#
# def reduction_D(D):
#     R = list(D)
#     list_low = [low_j(R, j) for j in range(len(R))]
#     for j in range(len(R)):
#         ilexiste = True
#         j0 = 0
#         while not (ilexiste) and j0 < j:
#             ilexiste = False
#             if list_low[j] == list_low[j0]:
#                 R[j] = sum_mod_2(R, j, j0)
#                 ilexiste = True
#             else:
#                 j0 += 1
#     return (R)
#
#
# def paires_persistence(R, lzyCplx, nbr_tem):
#     list_low = [low_j(R, j) for j in range(len(R))]
#     paires = []
#     for i in range(len(list_low)):
#         if list_low[i] != 0:
#             data, time, type = lzyCplx[i - nbr_tem]
#             paires.append([(list_low[i], i), type])
#     return (paires)
