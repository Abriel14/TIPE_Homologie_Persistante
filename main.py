import numpy as np


class Rips_complex:
    def __init__(self, points):
        """
        :type points: list
        """
        self.splxs = []
        self.nbr_splxs = len(self.splxs)
        self.zero_splxs = np.array(points)
        self.one_splxs = np.array([])
        self.two_splxs = np.array([])
        self.nbr_0_splxs = len(self.zero_splxs)
        self.nbr_1_splxs = len(self.one_splxs)
        self.nbr_2_splxs = len(self.two_splxs)
        self.D = self.max_distance(points)
        self.step = self.min_distance(points) / 2
        self.dates = [(0, 0)] * len(points)
        self.neighbours_matrix = np.array([])
        self.pers_diag = self.persistant_homology(self.zero_splxs, self.D, self.step)

    def max_distance(self, points):
        dist_max = 0
        for (a, b) in (points, points):
            dist = self.distance(a, b)
            if dist_max < dist:
                dist_max = dist
        return dist_max

    def distance(self, a, b):
        xa, ya = a
        xb, yb = b
        return (np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb)))

    def min_distance(self, points):
        dist_min = 0
        for (a, b) in (points, points):
            dist = self.distance(a, b)
            if dist_min > dist and a != b:
                dist_min = dist
        return dist_min

    def construct_network(self, dist_max, step):
        r = 0
        n = self.nbr_0_splxs
        for k in range(n):
            self.splxs.append((0, k))
        while r <= dist_max:
            for i in range(n):
                for j in range(i + 1, n):
                    a = self.zero_splxs[i]
                    b = self.zero_splxs[j]
                    if self.distance(a, b) < r:
                        np.append(self.one_splxs, ((i, j), self.nbr_splxs))
                        self.splxs.append((1, self.nbr_1_splxs))
                        self.dates.append((r, r))
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
                                        self.dates.append((r, r))
                                else:
                                    if e == d and f == b:
                                        np.append(self.two_splxs, (i, j, k))
                                        self.splxs.append((2, self.nbr_2_splxs))
                                        self.dates.append((r, r))
        return ()

    def construct_neighbours_matrix(self):
        n = self.nbr_splxs
        self.neighbours_matrix = np.zeros((n, n))
        for k in range(n):
            type, simplex = self.splxs[k]
            if type == 1:
                i, j = self.one_splxs[simplex]
                self.neighbours_matrix[i, k] = 1
                self.neighbours_matrix[j, k] = 1
            if type == 2:
                i, j, k = self.two_splxs[simplex]
                self.neighbours_matrix[self.one_splxs[i][1], k] = 1
                self.neighbours_matrix[self.one_splxs[j][1], k] = 1
                self.neighbours_matrix[self.one_splxs[k][1], k] = 1
        return ()

    def persistant_homology(self, points, dist_max: object, step: object) -> object:

        return ()
