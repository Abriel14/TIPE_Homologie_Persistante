import numpy as np


class Rips_complex:
    def __init__(self, points):
        """
        :type points: list
        """
        self.splxs = []
        self.points = points
        self.nbr_splxs = len(self.splxs)
        self.zero_splxs = np.array(points)
        self.one_splxs = np.array([])
        self.two_splxs = np.array([])
        self.nbr_0_splxs = len(self.zero_splxs)
        self.nbr_1_splxs = len(self.one_splxs)
        self.nbr_2_splxs = len(self.two_splxs)
        self.D = self.max_distance(points)
        self.step = self.min_distance(points) / 2
        self.birth_dates = [0] * len(points)
        self.connected_components_dates = []
        self.cycles_dates = []
        self.holes_dates = []
        self.neighbours_matrix = np.array([])
        self.homology_matrix = np.array([])
        self.pers_diag = self.execute_homology

    def distance(self, a, b):
        """        
        :param a: int tuple
        :param b: int tuple
        :return: eulerian distance between the two points
        """
        xa, ya = a
        xb, yb = b
        return (np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)))

    def max_distance(self, points):
        """
        :param points: np.array()
        :return: maximum distance between 2 points in a scatter plots
        """
        dist_max = 0
        for (a, b) in (points, points):
            dist = self.distance(a, b)
            if dist_max < dist:
                dist_max = dist
        return dist_max

    def min_distance(self, points):
        """
        :param points: np.array()
        :return: minimum distance between 2 points in a scatter plots
        """
        dist_min = 0
        for (a, b) in (points, points):
            dist = self.distance(a, b)
            if dist_min > dist and a != b:
                dist_min = dist
        return dist_min

    def construct_network(self):
        """
        construct the Ripps network for each distance r from 0 to dist_max with a given step
        :param dist_max: 
        :param step: 
        :return: void
        """
        r = 0
        n = self.nbr_0_splxs
        for k in range(n):
            self.splxs.append((0, k))
        while r <= self.D:
            # Create the 1_splxs
            for i in range(n):
                for j in range(i + 1, n):
                    a = self.zero_splxs[i]
                    b = self.zero_splxs[j]
                    if self.distance(a, b) < r:
                        np.append(self.one_splxs, ((i, j), self.nbr_splxs))
                        self.splxs.append((1, self.nbr_1_splxs))
                        self.birth_dates.append(r)
            # Create the 2_splxs
            for i in range(self.nbr_1_splxs):
                a, b = self.one_splxs[i][0]
                for j in range(self.nbr_1_splxs):
                    if i != j:
                        c, d = self.one_splxs[j][0]
                        if c == a:
                            for k in range(self.nbr_1_splxs):
                                e, f = self.one_splxs[k][0]
                                if b < d:
                                    if e == b and f == d:
                                        np.append(self.two_splxs, (i, j, k))
                                        self.splxs.append((2, self.nbr_2_splxs))
                                        self.birth_dates.append(r)
                                else:
                                    if e == d and f == b:
                                        np.append(self.two_splxs, (i, j, k))
                                        self.splxs.append((2, self.nbr_2_splxs))
                                        self.birth_dates.append(r)
                r += self.step
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
            type, simplex = self.splxs[k]
            # two possible conditions depending on the type of the simplex :
            # a 1_splx has 2 got edges whereas a 2_splx has got 3 of them.
            if type == 1:
                i, j = self.one_splxs[simplex]
                self.neighbours_matrix[i, k] = 1
                self.neighbours_matrix[j, k] = 1
            if type == 2:
                i, j, l = self.two_splxs[simplex]
                self.neighbours_matrix[self.one_splxs[i][1], k] = 1
                self.neighbours_matrix[self.one_splxs[j][1], k] = 1
                self.neighbours_matrix[self.one_splxs[l][1], k] = 1
        return ()

    @property
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
                    if low(j0, self.homology_matrix) == low(j, self.homology_matrix):
                        self.homology_matrix[:j, j] = (self.homology_matrix[:j, j0] + self.homology_matrix[:j, j]) % 2
                        test = True

        # Calculate the date of death of the connected components
        # for 0_simplexes, cycles for 1_splx ans holes for 2-splx
        for j in range(n):
            low_j = low(j, self.homology_matrix)
            if low_j >= 0 and self.homology_matrix[:j, low_j] == np.zeros((1, j)):
                birth = self.birth_dates[low_j]
                death = self.birth_dates[j]
                type, simplex = self.splxs[low_j]
                if type == 0:
                    self.connected_components_dates.append((birth, death))
                if type == 1:
                    self.cycles_dates.append((birth, death))
                if type == 2:
                    self.holes_dates.append((birth, death))
        return ()

    def execute_homology(self):
        self.construct_network()
        self.construct_neighbours_matrix()
        self.persistent_homology()
        return()

