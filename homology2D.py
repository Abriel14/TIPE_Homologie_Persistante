import numpy as np
from matplotlib.pyplot import *
import math as m
from operator import itemgetter
import random as rd
import math


class Point:
    """ Point class represents and manipulates x,y coords. """

    def __init__(self, x, y, numero):
        """ Create a new point en x,y """
        self.x = x
        self.y = y
        self.numero = numero

    def dist(self, p):
        """calcule la distance entre un point et un autre"""
        d = m.sqrt((self.x - p.x) ** 2 + (self.y - p.y) ** 2)
        return d


def affiche(tabP, name, scatt=20):
    """affiche un nuage de point"""
    px = [tabP[i].x for i in range(len(tabP))]
    py = [tabP[i].y for i in range(len(tabP))]
    scatter(px, py, scatt)
    xlabel('x')
    ylabel('y')
    savefig(name, format='pdf', transparent=True)
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
    # simplexes est un grand tableau regroupant les informations [simplexe, date, type de simplexe] sur les simplexes créés
    simplexes.sort(key=itemgetter(1))
    # on le trie par date d'apparition des simplexes croissante
    print("Witness complex créé")
    return (simplexes, nbrL)


def calcule_D(lazycplx, nbrL):
    """calcule la matrice D"""
    D = np.zeros((len(lazycplx) + nbrL, len(lazycplx) + nbrL))
    pos_unsimplexe = np.zeros((nbrL, nbrL))
    list_low = np.zeros(len(lazycplx) + nbrL)
    for v in range(len(lazycplx)):
        pos = v + nbrL
        (splx, time, typpe) = lazycplx[v]
        if typpe == 1:
            # c'est un 1 simplexe
            pos_unsimplexe[splx[0], splx[1]] = pos
            pos_unsimplexe[splx[1], splx[0]] = pos
            D[splx[0], pos] = 1
            D[splx[1], pos] = 1
            list_low[pos] = max(splx[0], splx[1])

        else:
            # c'est un 2 simplexe
            D[int(pos_unsimplexe[splx[0], splx[1]]), pos] = 1
            D[int(pos_unsimplexe[splx[1], splx[2]]), pos] = 1
            D[int(pos_unsimplexe[splx[0], splx[2]]), pos] = 1
            list_low[pos] = max(int(pos_unsimplexe[splx[0], splx[1]]), int(pos_unsimplexe[splx[1], splx[2]]),
                                int(pos_unsimplexe[splx[0], splx[2]]))
    print("Matrice de l'application bord initialisée")
    return (D, list_low)


def low(mat, j):
    """calcule low(j), low(j) est le plus grand indice i tel que la i-ème valaur de la colonne j de mat vaut 1 """
    res = 0
    for i in range(len(mat)):
        if mat[i][j] == 1:
            res = i
    return (res)


def reduction_D(mat,listlow):
    """procède à la réduction de la matrice D"""
    for j in range(len(mat)):
        ilexiste = True
        while listlow[j] != 0 and ilexiste:
            ilexiste = False
            j0 = 0
            while j0 < j and not (ilexiste):
                if listlow[j0] == listlow[j]:
                    ilexiste = True
                    mat[:, j] = (mat[:, j] + mat[:, j0]) % 2
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
            res.append((sorted([listlow[i], i]), type))
    return (res)


def diag_pers(paires, name):
    h0_birth = []
    h0_death = []
    h1_birth = []
    h1_death = []
    for p in paires:
        [[b, d], type] = p
        if type == 1:
            h0_birth.append(b)
            h0_death.append(d)
        if type == 2:
            h1_birth.append(b)
            h1_death.append(d)
    figure(figsize=(10, 10), dpi=90)
    subplot(211)
    title("H0")
    xlabel('Rayon epsilon')
    ylabel("Indice de la paire")
    hlines(range(len(h0_birth)), h0_birth, h0_death, colors='b')
    subplot(212)
    title("H1")
    xlabel('Rayon epsilon')
    ylabel("Indice de la paire")
    hlines(range(len(h1_birth)), h1_birth, h1_death, colors='b')
    savefig(name, format='pdf', transparent=True)
    show()

def hom_pers(P,witness_pourc,file_name, v = 1):
    cplx, nbL = Lazy_Witness_Complex(P, witness_pourc, v)
    D, list_low = calcule_D(cplx, nbL)
    DD = reduction_D(D, list_low)
    C = paires_pers(D, cplx, nbL)
    diag_pers(C, file_name)


def annulus(N, r1, r2):
    anneau = []
    k = 0
    while k< N:
        p = Point(-1 + 2 * rd.random(), -1 + 2 * rd.random(), k)
        d = math.sqrt((p.x) ** 2 + (p.y) ** 2)
        if d > r1 and d < r2:
            anneau.append(p)
            k += 1
    return anneau


def limite(Nmin, Nmax, r1, r2, nbr_test):
    res = []
    for N in range(Nmax, Nmin, -1):
        print(N)
        general_mean = 0
        for i in range(nbr_test):
            anneau = annulus(N, r1, r2)
            cplx, nbL = Lazy_Witness_Complex(anneau, 0.1, 1)
            D = calcule_D(cplx, nbL)
            DD = reduction_D(D)
            C = paires_pers(D, cplx, nbL)
            k0 = 0
            highest = 0
            for k in range(len(C)):
                (time_b, time_d), type = C[k]
                length = time_d - time_b
                if length > highest:
                    highest = length
                    k0 = k
            mean_low = 0
            for k in range(len(C)):
                (time_b, time_d), type = C[k]
                if k != k0:
                    mean_low += time_d - time_b
            mean_low = mean_low / (len(C) - 1)
            general_mean += abs(highest / mean_low)
        general_mean = general_mean / nbr_test
        res.append(general_mean)
    plot(range(Nmax, Nmin, -1), res)
    show()
