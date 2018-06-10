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


def affiche(tabP, name, scatt=20):
    """affiche un nuage de point"""
    px = [tabP[i].x for i in range(len(tabP))]
    py = [tabP[i].y for i in range(len(tabP))]
    pz = [tabP[i].z for i in range(len(tabP))]
    fig = figure(figsize=(7, 6.5))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)
    ax.scatter(px, py, pz)
    savefig(name, format='pdf', transparent=True)
    show()


def mini_maxi(tabP, delta, i, j, nbrL):
    """calcule le minimum des max(d(a1,w),d(a2,w)) - delta[w] pour tous les point w n'appartenants pas à la filtration"""
    res = [abs(max(tabP[w].dist(tabP[i]), tabP[w].dist(tabP[j])) - delta[w - nbrL]) for w in range(nbrL, len(tabP))]
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
    print("Witness complex créé")
    return (simplexes, nbrL)


def calcule_D(lazycplx, nbrL):
    """calcule la matrice D"""
    D = [[] for i in range(len(lazycplx) + nbrL)]
    pos_unsimplexe = np.zeros((nbrL, nbrL))
    pos_deuxsimplexe = np.zeros((nbrL, nbrL, nbrL))
    list_low = np.zeros(len(lazycplx) + nbrL)
    for v in range(len(lazycplx)):
        pos = v + nbrL
        (splx, time, type) = lazycplx[v]
        if type == 1:
            # c'est un 1 simplexe
            pos_unsimplexe[splx[0], splx[1]] = pos
            pos_unsimplexe[splx[1], splx[0]] = pos
            D[pos] = sorted([splx[0], splx[1]])
            # D[splx[0], pos] = 1
            # D[splx[1], pos] = 1
            list_low[pos] = max(splx[0], splx[1])
        elif type == 2:
            # c'est un 2 simplexe
            pos_deuxsimplexe[splx[0], splx[1], splx[2]] = pos
            pos_deuxsimplexe[splx[0], splx[2], splx[1]] = pos
            pos_deuxsimplexe[splx[1], splx[0], splx[2]] = pos
            pos_deuxsimplexe[splx[1], splx[2], splx[0]] = pos
            pos_deuxsimplexe[splx[2], splx[1], splx[0]] = pos
            pos_deuxsimplexe[splx[2], splx[0], splx[1]] = pos
            D[pos] = sorted([int(pos_unsimplexe[splx[0], splx[1]]), int(pos_unsimplexe[splx[1], splx[2]]),
                             int(pos_unsimplexe[splx[0], splx[2]])])
            # D[int(pos_unsimplexe[splx[0], splx[1]]), pos] = 1
            # D[int(pos_unsimplexe[splx[1], splx[2]]), pos] = 1
            # D[int(pos_unsimplexe[splx[0], splx[2]]), pos] = 1
            list_low[pos] = max(int(pos_unsimplexe[splx[0], splx[1]]), int(pos_unsimplexe[splx[1], splx[2]]),
                                int(pos_unsimplexe[splx[0], splx[2]]))
        else:
            # c'est un 3 simplexe
            D[pos] = sorted(
                [int(pos_deuxsimplexe[splx[0], splx[1], splx[2]]), int(pos_deuxsimplexe[splx[0], splx[1], splx[3]]),
                 int(pos_deuxsimplexe[splx[0], splx[2], splx[3]]), int(pos_deuxsimplexe[splx[1], splx[2], splx[3]])])
            # D[int(pos_deuxsimplexe[splx[0], splx[1], splx[2]]), pos] = 1
            # D[int(pos_deuxsimplexe[splx[0], splx[1], splx[3]]), pos] = 1
            # D[int(pos_deuxsimplexe[splx[0], splx[2], splx[3]]), pos] = 1
            # D[int(pos_deuxsimplexe[splx[1], splx[2], splx[3]]), pos] = 1
            list_low[pos] = max(int(pos_deuxsimplexe[splx[0], splx[1], splx[2]]),
                                int(pos_deuxsimplexe[splx[0], splx[1], splx[3]]),
                                int(pos_deuxsimplexe[splx[0], splx[2], splx[3]]),
                                int(pos_deuxsimplexe[splx[1], splx[2], splx[3]]))
    print("Matrice de l'application bord initialisée")
    return (D, list_low)


def low(mat, j):
    """calcule low(j), low(j) est le plus grand indice i tel que la i-ème valaur de la colonne j de mat vaut 1 """
    if mat[j] == []:
        return 0
    else:
        return mat[j][-1]

def fusion(T1,T2) :
    if T1==[] :return T2
    if T2==[] :return T1
    if T1[0]<T2[0] :
        return [T1[0]]+fusion(T1[1 :],T2)
    elif T1[0]==T2[0]:
        return fusion(T1[1:], T2[1:])
    else:
        return [T2[0]] + fusion(T1, T2[1:])


def reduction_D(mat, listlow):
    """procède à la réduction de la matrice D"""
    # listlow = [low(mat, j) for j in range(len(mat))]
    for j in range(len(mat)):
        print(j / len(mat))
        ilexiste = True
        while listlow[j] != 0 and ilexiste:
            ilexiste = False
            j0 = 0
            while j0 < j and not (ilexiste):
                if listlow[j0] == listlow[j]:
                    ilexiste = True
                    mat[j] = fusion(mat[j], mat[j0])
                    listlow[j] = low(mat, j)
                else:
                    j0 = j0 + 1
    print("Matrice de l'application bord réduite")
    return (mat)


def paires_pers(D, cplx, nbrL):
    listlow = [low(D, j) for j in range(len(D))]
    res = []
    for i in range(len(listlow)):
        if listlow[i] != 0:
            (splx, time_d, type) = cplx[i - nbrL]
            (splx, time_b, typeb) = cplx[listlow[i] - nbrL]
            res.append(((listlow[i], i), type))
    return (res)


def diag_pers(paires, name):
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
    figure(figsize=(10, 15), dpi=90)
    subplot(311)
    title("H0")
    xlabel('Rayon epsilon')
    ylabel("Indice de la paire")
    hlines(range(len(h0_birth)), h0_birth, h0_death, colors='b')
    subplot(312)
    title("H1")
    xlabel('Rayon epsilon')
    ylabel("Indice de la paire")
    hlines(range(len(h1_birth)), h1_birth, h1_death, colors='b')
    subplot(313)
    title("H2")
    xlabel('Rayon epsilon')
    ylabel("Indice de la paire")
    hlines(range(len(h2_birth)), h2_birth, h2_death, colors='b')
    savefig(name, format='pdf', transparent=True)
    show()


def hom_pers(P, witness_pourc, file_name, v=1):
    cplx, nbL = Lazy_Witness_Complex(P, witness_pourc, v)
    D, list_low = calcule_D(cplx, nbL)
    DD = reduction_D(D, list_low)
    C = paires_pers(D, cplx, nbL)
    diag_pers(C, file_name)
