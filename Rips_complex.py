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
        self.one_splxs = []  # list of all 1-simplexes
        self.two_splxs = []  # list of all 2-simplexes
        self.nbr_0_splxs = len(self.zero_splxs)
        self.nbr_1_splxs = len(self.one_splxs)
        self.nbr_2_splxs = len(self.two_splxs)
        self.D = self.max_distance()  # The maximum distance in the scatter plot
        self.step = self.min_distance() / 3  # the step used in the construction of the simplex
        self.birth_dates = [0] * len(self.points)  # list of the date of creation of each simplexes
        self.connected_components_dates = []  # list of the dates of 0-simplexes's invariants ruptures
        self.cycles_dates = []  # list of the dates of 1-simplexes's invariants ruptures
        self.holes_dates = []  # list of the dates of 2-simplexes's invariants ruptures
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

    def construct_network(self):
        """
        construct the Ripps network for each distance r from 0 to dist_max with a given step
        :return: void
        """
        r = 0
        n = self.nbr_0_splxs
        for k in range(n):
            self.splxs.append((0, (0, k)))
            self.nbr_splxs+=1
        while r <= (self.D + self.step):
            # Create the 1_splxs
            for i in range(n):
                for j in range(i + 1, n):
                    a = self.zero_splxs[i]
                    b = self.zero_splxs[j]
                    if r - self.step <= self.distance(a, b) <= r:
                        self.one_splxs.append(((i, j),self.nbr_splxs))
                        self.splxs.append((1, self.nbr_1_splxs))
                        self.birth_dates.append(r)
                        self.nbr_1_splxs+=1
                        self.nbr_splxs+=1
            # Create the 2_splxs
            if self.nbr_1_splxs != 0:
                for i in range(self.nbr_1_splxs):
                    a, b = self.one_splxs[i][0]
                    for j in range(i,self.nbr_1_splxs):
                        if i != j:
                            c, d = self.one_splxs[j][0]
                            if c == a:
                                for k in range(self.nbr_1_splxs):
                                    e, f = self.one_splxs[k][0]
                                    if max((self.distance((a,b),(c,d)),self.distance((a,b),(e,f)),self.distance((c,d),(e,f))))>= r - self.step:
                                        if b < d:
                                            if e == b and f == d:
                                                self.two_splxs.append((i, j, k))
                                                self.splxs.append((2, self.nbr_2_splxs))
                                                self.birth_dates.append(r)
                                                self.nbr_2_splxs +=1
                                                self.nbr_splxs += 1
                                        else:
                                            if e == d and f == b:
                                                self.two_splxs.append((i, j, k))
                                                self.splxs.append((2, self.nbr_2_splxs))
                                                self.birth_dates.append(r)
                                                self.nbr_2_splxs +=1
                                                self.nbr_splxs += 1
            r += self.step
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
            low_j = -1
            for k in range(j):
                if R[k, j] == 1:
                    low_j = k
            return (low_j)

        self.homology_matrix = self.neighbours_matrix.copy()
        n = self.nbr_splxs

        # Apply the persistence algorithm

        for j in range(1, n):
            test = True
            while test == True:
                test = False
                for j0 in range(j):
                    if low(j0, self.homology_matrix) == low(j, self.homology_matrix) and low(j0, self.homology_matrix) !=-1:
                        self.homology_matrix[:j, j] = (self.homology_matrix[:j, j0] + self.homology_matrix[:j, j]) % 2
                        test = True
                        # if all(self.homology_matrix[:j, j]) == 0:
                        #     type, simplex = self.splxs[j]
                        #     birth = self.birth_dates[j0]
                        #     death = self.birth_dates[j]
                        #     if type == 0:
                        #         self.connected_components_dates.append((birth, death))
                        #     if type == 1:
                        #         self.cycles_dates.append((birth, death))
                        #     if type == 2:
                        #         self.holes_dates.append((birth, death))
                        # print(self.homology_matrix)
        # Calculate the date of death of the connected components
        # for 0_simplexes, cycles for 1_splx ans holes for 2-splx
        for j in range(n):
            low_j = low(j, self.homology_matrix)
            if low_j >= 0 and self.homology_matrix[j:, low_j].all() == 0:
                birth = self.birth_dates[low_j]
                death = self.birth_dates[j]
                type, simplex = self.splxs[j]
                if type == 0:
                    self.connected_components_dates.append((birth, death))
                if type == 1:
                    self.cycles_dates.append((birth, death))
                if type == 2:
                    self.holes_dates.append((birth, death))
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
        print(self.connected_components_dates)
        print(self.cycles_dates)
        print(self.holes_dates)
        # figure()
        # hlines(range(len(self.connected_components_dates)),self.connected_components_dates[:][0], self.connected_components_dates[:][1] )
        # title('Connected components')
        figure()
        hlines(range(len(self.cycles_dates)),self.cycles_dates[:][0], self.cycles_dates[:][1])
        title('Cycles')
        figure()
        hlines(range(len(self.holes_dates)), self.holes_dates[:][0],self.holes_dates[:][1])
        title('holes')
        show()
        return ()
