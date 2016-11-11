import numpy as np
import paraFEM

u_knots = np.array([ 0., 0., 1., 1.])
v_knots = np.array([ 0., 0., 1., 1.])
u_degreee = 1
v_degree = 1
weights = np.array([1., 1., 1., 1.])

a = paraFEM.NurbsBase(u_knots, v_knots, weights, u_degreee, v_degree)
a.computeFirstDerivatives()
u = np.linspace(u_knots[0], u_knots[-1], 21)
v = np.linspace(v_knots[0], v_knots[-1], 21)

uv = [[i, j] for i in u for j in v]
points = a.getDuMatrix(np.array(uv)) * np.array([0, 1, 1, 0])
print(points)