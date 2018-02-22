import numpy as np
from matplotlib.pyplot import *


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
        # self.birth_dates = [0] * len(self.points)  # list of the date of creation of each simplexes
        self.connected_components_birth = []  # list of the birth dates of 0-simplexes's invariants ruptures
        self.connected_components_death = []  # list of the death dates of 0-simplexes's invariants ruptures
        self.pers_pairs_birth = []  # list of the biths of persistent pairs
        self.pers_pairs_death = []  # list of the death of persistent pairs
        self.low_j_to_j_list = []
        self.neighbours_matrix = np.array([])
        self.homology_matrix = np.array([])
        self.pers_diag = self.execute_homology()

    def distance(self, a, b):
        """        
        :param a: int tuple
        :param b: int tuple
        :return: eulerian distance between the two points
        """
        xa, ya = a
        xb, yb = b
        return (np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)))

    def max_distance(self):
        """
        :param points: np.array()
        :return: maximum distance between 2 points in a scatter plots
        """
        dist_max = 0
        for a in self.points:
            for b in self.points:
                dist = self.distance(a, b)
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
                dist = self.distance(a, b)
                if dist_min > dist and a != b:
                    dist_min = dist
        return dist_min

    def find_next_edge(self, last_dist):
        dist_min = self.D
        edge = -1, -1
        for i in range(self.nbr_0_splxs):
            for j in range(i, self.nbr_0_splxs):
                a = self.points[i]
                b = self.points[j]
                dist = self.distance(a, b)
                if dist_min >= dist > last_dist:
                    dist_min = dist
                    edge = (i, j)
        return (dist_min, edge)

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
            # find the next edge to treat
            r, edge = self.find_next_edge(r)
        print("Network created")
        return ()

    def construct_neighbours_matrix(self):
        """
        Constructs the neighbours matrix:
        neighbour_matrix[i,k] is equal to 1 iif the simplex indexed by k in
        the splxs matrix is an edge of the simplex indexed by i
        :return: void
        """
        n = self.nbr_splxs
        self.neighbours_matrix = np.zeros((n, n))
        for k in range(n):
            type, simplex_index = self.splxs[k]
            # two possible conditions depending on the type of the simplex :
            # a 1_splx has 2 got edges whereas a 2_splx has got 3 of them.
            if type == 1:
                i, j = self.one_splxs[simplex_index][0]
                self.neighbours_matrix[i, k] = 1
                self.neighbours_matrix[j, k] = 1
            if type == 2:
                i, j, l = self.two_splxs[simplex_index]
                self.neighbours_matrix[self.one_splxs[i][1], k] = 1
                self.neighbours_matrix[self.one_splxs[j][1], k] = 1
                self.neighbours_matrix[self.one_splxs[l][1], k] = 1

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
            low_j = 0
            for k in range(j):
                if R[k, j] == 1:
                    low_j = k
            return (low_j)

        self.homology_matrix = self.neighbours_matrix.copy()
        n = self.nbr_0_splxs
        N = self.nbr_splxs
        # initilize the low_j matrix
        self.low_j_to_j_list = N * [0]
        # Apply the persistence algorithm
        j = 0
        while low(j, self.homology_matrix) == 0:
            j+=1
        self.low_j_to_j_list[low(j, self.homology_matrix)] = j
        j+=1
        while j<N:
            low_j = low(j,self.homology_matrix)
            j0 = self.low_j_to_j_list[low_j]
            while j0 != 0:
                self.homology_matrix[:j, j] = (self.homology_matrix[:j, j0] + self.homology_matrix[:j, j]) % 2
                low_j = low(j, self.homology_matrix)
                j0 = self.low_j_to_j_list[low_j]
            if low_j !=0:
                self.low_j_to_j_list[low_j] = j
            j+=1
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
            if low_j != 0:
                # print(low_j,j)
                self.pers_pairs_birth.append(low_j)
                self.pers_pairs_death.append(j)

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
        hlines(range(0, -len(self.pers_pairs_death), -1), self.pers_pairs_birth, self.pers_pairs_death)
        show()
        return ()
