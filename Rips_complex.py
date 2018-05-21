import numpy as np
from matplotlib.pyplot import *


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

def affiche(tabP,scatt = 20):
    px = [tabP[i].x for i in range(len(tabP))]
    py = [tabP[i].y for i in range(len(tabP))]
    scatter(px,px,scatt)
    xlabel('x')
    ylabel('y')
    show()



class Rips_complex:
    def __init__(self, points):
        """
        :type points: list
        """
        self.splxs = []  # full list of all the simplexes
        self.points = points  # list of all points of the simplex
        self.nbr_splxs = len(self.splxs)
        self.zero_splxs = self.points  # list of all 0-simplexes
        self.one_splxs = []  # list of all 1-simplexes represented by a pair of 0-splx
        self.two_splxs = []  # list of all 2-simplexes represented by a triplet of 1-splx
        self.nbr_0_splxs = len(self.zero_splxs)
        self.nbr_1_splxs = len(self.one_splxs)
        self.nbr_2_splxs = len(self.two_splxs)
        self.D = self.max_distance()  # The maximum distance in the scatter plot
        self.d = self.min_distance()
        self.dist_appearance = [0]*self.nbr_0_splxs
        # self.birth_dates = [0] * len(self.points)  # list of the date of creation of each simplexes
        self.connected_components_birth = []  # list of the birth dates of 0-simplexes's invariants ruptures
        self.connected_components_death = []  # list of the death dates of 0-simplexes's invariants ruptures
        self.pers_pairs_birth = []  # list of the biths of persistent pairs
        self.pers_pairs_death = []  # list of the death of persistent pairs
        self.h0_birth = []
        self.h0_death = []
        self.h1_birth = []
        self.h1_death = []
        self.low_j_to_j_list = []
        self.neighbours_matrix = []
        self.homology_matrix = []
        self.pers_diag = self.execute_homology()



    def max_distance(self):
        """
        :param points: np.array()
        :return: maximum distance between 2 points in a scatter plots
        """
        dist_max = 0
        for a in self.points:
            for b in self.points:
                dist = a.dist(b)
                if dist_max < dist:
                    dist_max = dist
        return dist_max

    def min_distance(self):
        """
        :param points: np.array()
        :return: minimum distance between 2 points in a scatter plots
        """
        dist_min = self.D
        for a in self.points:
            for b in self.points:
                dist = a.dist(b)
                if dist_min > dist and a != b:
                    dist_min = dist
        return dist_min

    def mean_dist(self):
        dist_sum = 0
        for a in self.points:
            for b in self.points:
                if a != b:
                    dist_sum += a.dist(b)
        return (dist_sum / (len(self.points) * (len(self.points) - 1)))

    def find_next_edge(self, last_dist):
        dist_min = self.D
        edge = -1, -1
        for i in range(self.nbr_0_splxs):
            for j in range(i, self.nbr_0_splxs):
                a = self.points[i]
                b = self.points[j]
                dist = a.dist(b)
                if dist_min >= dist > last_dist:
                    dist_min = dist
                    edge = (i, j)
        # if dist_min > self.D / 8.0:
        #     edge = -1, -1
        return (dist_min, edge)

    def sum_column(self, i, j, M):
        column_i = sorted(M[i])
        column_j = sorted(M[j])
        if len(column_j) < len(column_i):
            temp = column_j
            column_j = column_i
            column_i = temp
        sum = []
        while column_j != []:
            if column_i != []:
                a = column_i[-1]
                b = column_j[-1]
                if a < b:
                    sum.append(column_j.pop())
                if a > b:
                    sum.append(column_i.pop())
                if a == b:
                    column_i.pop()
                    column_j.pop()
            else:
                sum.append(column_j.pop())
        return (sorted(sum))

    def construct_network(self):
        """
        construct the Ripps network for each distance r from 0 to dist_max with a given step
        :return: void
        """
        r = 0
        n = self.nbr_0_splxs
        for k in range(n):
            self.splxs.append((0, (0, k)))
            self.nbr_splxs += 1
        r, edge = self.find_next_edge(r)
        # this while loop finds the new edge to treat and add it to the 1-splx list and then finds out if a 2-splx is created
        while edge != (-1, -1):
            # Add the new edge
            self.one_splxs.append((edge, self.nbr_splxs))
            self.splxs.append((1, self.nbr_1_splxs))
            self.nbr_1_splxs += 1
            self.nbr_splxs += 1
            self.dist_appearance.append(r)
            a, b = edge
            # find out if a 2-splx has been created
            for i in range(self.nbr_1_splxs - 1):
                c, d = self.one_splxs[i][0]
                if d == a:
                    for j in range(i + 1, self.nbr_1_splxs - 1):
                        e, f = self.one_splxs[j][0]
                        if e == c and f == b:
                            self.two_splxs.append((self.nbr_1_splxs - 1, i, j))
                            self.splxs.append((2, self.nbr_2_splxs))
                            self.nbr_2_splxs += 1
                            self.nbr_splxs += 1
                            self.dist_appearance.append(r)
            # find the next edge to treat
            r, edge = self.find_next_edge(r)
        print("Network created")
        return ()

    def construct_neighbours_matrix(self):
        """
        Constructs the neighbours matrix:
        neighbour_matrix[i,k] is equal to 1 iif the simplex indexed by i in
        the splxs matrix is an edge of the simplex indexed by k
        :return: void
        """
        N = self.nbr_splxs
        self.neighbours_matrix = N * [[]]
        for k in range(N):
            type, simplex_index = self.splxs[k]
            # two possible conditions depending on the type of the simplex :
            # a 1_splx has 2 edges whereas a 2_splx has 3.
            if type == 1:
                i, j = self.one_splxs[simplex_index][0]
                self.neighbours_matrix[k] = sorted([i, j])
                # self.neighbours_matrix[k].append(i)
                # self.neighbours_matrix[k].append(j)
            if type == 2:
                i, j, l = self.two_splxs[simplex_index]
                self.neighbours_matrix[k] = sorted([self.one_splxs[i][1], self.one_splxs[j][1], self.one_splxs[l][1]])
                # self.neighbours_matrix[self.one_splxs[i][1], k] = 1
                # self.neighbours_matrix[self.one_splxs[j][1], k] = 1
                # self.neighbours_matrix[self.one_splxs[l][1], k] = 1
        print("Neighbours matrix created")
        return ()

    def persistent_homology(self):
        """
        function that calculates the persistent homology of the scatter plot
        and adds the dates of birth and death to the k_splx_dates list
        depending on the type of the simplex k that is filled with an edge  
        :return: void
        """

        def low(j, R):
            """
            :return: maximum line index of the column j in the matrix R with a 1 in it
            """
            if R[j] == []:
                return (-1)
            else:
                return (sorted(R[j])[-1])

                # low_j = 0
                # for k in range(j):
                #     if R[k, j] == 1:
                #         low_j = k
                # return (low_j)

        N = self.nbr_splxs
        self.homology_matrix = self.neighbours_matrix[:]
        n = self.nbr_0_splxs
        # initilize the low_j matrix
        self.low_j_to_j_list = N * [-1]
        # Apply the persistence algorithm
        j = 0
        while low(j, self.homology_matrix) == -1:
            j += 1
        self.low_j_to_j_list[low(j, self.homology_matrix)] = j
        j += 1
        while j < N:
            low_j = low(j, self.homology_matrix)
            j0 = self.low_j_to_j_list[low_j]
            while j0 != -1:
                self.homology_matrix[j] = self.sum_column(j, j0, self.homology_matrix)
                # self.homology_matrix[:j, j] = (self.homology_matrix[:j, j0] + self.homology_matrix[:j, j]) % 2
                low_j = low(j, self.homology_matrix)
                j0 = self.low_j_to_j_list[low_j]
            if low_j != -1:
                self.low_j_to_j_list[low_j] = j
            j += 1
            if j % 10 == 0:
                print(j / N)
        # for j in range(1, N):
        #     test = True
        #     while test:
        #         test = False
        #         for j0 in range(j):
        #             if low(j0, self.homology_matrix) == low(j, self.homology_matrix) \
        #                     and low(j0, self.homology_matrix) != 0:
        #                 self.homology_matrix[:j, j] = (self.homology_matrix[:j, j0] + self.homology_matrix[:j, j]) % 2
        #                 test = True
        #     if j % 10 == 0:
        #         print(np.log(j + 1) / np.log(N))

        for j in range(N):
            low_j = low(j, self.homology_matrix)
            if low_j != -1:
                # print(low_j,j)
                # self.pers_pairs_birth.append(self.dist_appearance[low_j])
                # self.pers_pairs_death.append(self.dist_appearance[j])
                if self.splxs[low_j][0] == 0:
                    self.h0_birth.append(self.dist_appearance[low_j])
                    self.h0_death.append(self.dist_appearance[j])
                    print(low_j)
                else:
                    self.h1_birth.append(self.dist_appearance[low_j])
                    self.h1_death.append(self.dist_appearance[j])
        print("persistant homology achieved")
        return ()

    def execute_homology(self):
        self.construct_network()
        self.construct_neighbours_matrix()
        self.persistent_homology()
        return ()

    def show_figure(self):
        figure()
        plot(self.zero_splxs[:][0], self.zero_splxs[:][1])
        show()
        print(self.D)
        print(self.splxs)
        print(self.zero_splxs)
        print(self.one_splxs)
        print(self.two_splxs)
        return ()

    def show_pers_diagram(self):
        # figure()
        # hlines(range(len(self.connected_components_birth)), self.connected_components_birth,
        #        self.connected_components_death)
        # title('Connected components')
        # figure()
        # hlines(range(len(self.holes_birth)), self.holes_birth, self.holes_death)
        # title('holes')

        figure()
        title("H0")
        hlines(range(len(self.h0_birth)), self.h0_birth, self.h0_death)
        figure()
        title("H1")
        hlines(range(len(self.h1_birth)), self.h1_birth, self.h1_death)
        show()
        return ()
