import numpy as np
from hjälpfunktioner import *
V, K = generate_mesh_2d(1, 1)

print(assemble_vector_2d(poisson_rhs, V, K))

# Konvertera till vanlig matris och printa





    # Compute area of triangle
  #  a = X[:, 1] - X[:, 0]
  #  b = X[:, 2] - X[:, 0]
  #  area = 0.5 * abs(a[0] * b[1] - a[1] * b[0])
 
    # Compute quadrature points
  #  points = zeros((2, 1))
  #  points[:, 0] = (X[:, 0] + X[:, 1] + X[:, 2]) / 3.0

    # Compute quadrature weights
   # weights = array([area])

   # return points, weights