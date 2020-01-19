class Case(object):
    def __init__(self, elements, constraints):
        self.elements = elements
        self.constraints = constraints # same shape as points 0: fixed, 1: no constraint

    def run_case(self):
        pass


def Element(object):
    def __init__(self, points, nodes, coord_sys=None):
        self.points = points
        self.nodes = nodes
        self.coord_sys = coord_sys or self.compute_coord_sys() # 3 x 3 array

    def add_local_K(self, K_global):
        # 1 agl
        # 2 tgl
        # 3 xl = agl@xg
        # 4 dn
        # 5 J = dN @ xl
        # 6 Bp = (J^-1)@dN
        # 6 B
        # 7 det(J)
        # 8 K_local
        pass

    @property
    def global_points(self):
        return self.points[self.nodes]

    @property
    def global_indices(self):
        n = np.array(self.nodes) * 3
        return np.array([n + 0, n + 1, n + 2]).T.flatten().tolist()

    def insert_local_K(self, K_local, K_global):
        arr = self.global_indices
        for il, ig in enumerate(arr):
            for jl, jg in enumerate(arr):
                K.append([ig, jg, K_local[il, jl]])


class Truss(Element):
    def add_local_K(self, K_global):
        # 1 agl
        agl = self.coord_sys
        # 2 tgl
        tgl = agl[0]
        # 3 xl = agl@xg
        xg = self.global_points
        xl = tgl @ xg
        xl -= np.mean(xl)
        # 4 dn
        dN = np.array([-1. , 1.])
        # 5 J = dN @ xl
        J = dN @ xl
        # 6 Bp = (J^-1)@dN
        Bp = np.linalg.inv(J) @ dN
        # 6 B
        B = Bp
        # 7 det(J) = J (scalar)
        # 8 K_local
        K_local = tgl.T @ B.T @ C @ B @ tgl * J
        print(K_local)
        pass

    def compute_coord_sys(self):
        p1, p2 = self.global_points
        t = p2 - p1
        t /= np.linalg.norm(t)
        n1 = np.random.rand(3)
        n2 = np.cross(t, n1)
        n2 /= np.linalg.norm(n2)
        n1 = np.cross(n2, t)
        return np.array(t, n1, n2)
        

if __name__ == "__main__":
    points = np.array([[1, 1, 1],
                       [2, 1, 1]])
    nodes = [0, 1]
    truss_element = Truss(points, nodes)
    K_global = []
    truss_element.add_local_K(K_global)