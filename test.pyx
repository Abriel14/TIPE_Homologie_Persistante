# noinspection PyUnresolvedReferences
cimport cython
cimport numpy as np

ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE_f

# noinspection PyPep8Naming,PyUnresolvedReferences
@cython.boundscheck(False)  # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function

def low(np.ndarray[DTYPE_t, ndim=2] mat,DTYPE_t j):
    """calcule low(j), low(j) est le plus grand indice i tel que la i-ème valaur de la colonne j de mat vaut 1 """
    res = -1
    for i in range(len(mat)):
        if mat[i][j] == 1:
            res = i
    return (res)
def reduction_D(np.ndarray[DTYPE_f, ndim=2] mat):
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