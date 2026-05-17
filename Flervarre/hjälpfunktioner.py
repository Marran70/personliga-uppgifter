from numpy import array, zeros, dot, cross, shape
from numpy.linalg import inv
from pylab import *
#5.1
def generate_mesh_2d(nx, ny):

    hx = 1 / nx
    hy = 1 / ny

    V = []
    K = []

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            x = ix * hx
            y = iy * hy
            V.append((x, y))

    for iy in range(ny):
        for ix in range(nx):
            v0 = iy * (nx + 1) + ix
            v1 = v0 + 1
            v2 = v0 + (nx + 1)
            v3 = v1 + (nx + 1)
            K.append((v0, v1, v3))
            K.append((v0, v3, v2))

    return array(V), array(K)

#5.3
def generate_edges_2d(K):

    E = set()

    for t in K:
        edges = [(t[0], t[1]), (t[0], t[2]), (t[1], t[2])]
        for e in edges:
            E.add(tuple(sorted(e)))

    return list(E)


#5.5
def plot_mesh_2d(V, K):

    for t in enumerate(K):
        for j in range(3):
            for k in range(j, 3):
                x = [V[t[j]][0], V[t[k]][0]]
                y = [V[t[j]][1], V[t[k]][1]]
                plot(x, y)

    for v in enumerate(V):
        plot(v[0], v[1], "o")





#5.7
def evaluate_basis_2d(i, x, X):

    # Compute Jacobian (linear part) of affine map
    J = zeros((2, 2))
    J[:, 0] = X[:, 1] - X[:, 0]
    J[:, 1] = X[:, 2] - X[:, 0]

    # Compute inverse of Jacobian
    Jinv = inv(J)

    # Map point to reference element
    xhat = dot(Jinv, x - X[:, 0])

    # Evaluate basis function
    if i == 0:
        phi = 1 - xhat[0] - xhat[1]
        grad_phi = dot(Jinv.T, (-1, -1))
    elif i == 1:
        phi = xhat[0]
        grad_phi = dot(Jinv.T, (1, 0))
    else:
        phi = xhat[1]
        grad_phi = dot(Jinv.T, (0, 1))

    return phi, grad_phi


#5.9
def generate_quadrature_2d(X):

# Compute area of triangle
    a = X[:, 1] - X[:, 0]
    b = X[:, 2] - X[:, 0]
    area = 0.5 * abs(a[0] * b[1] - a[1] * b[0])

    # Compute quadrature points (midpoints of the edges)
    # p0: midpoint between vertex 0 and 1
    # p1: midpoint between vertex 1 and 2
    # p2: midpoint between vertex 2 and 0
    p0 = (X[:, 0] + X[:, 1]) / 2.0
    p1 = (X[:, 1] + X[:, 2]) / 2.0
    p2 = (X[:, 2] + X[:, 0]) / 2.0

    # Stack the points into a 2x3 matrix
    points = np.column_stack((p0, p1, p2))

    # Compute quadrature weights
    # For this rule, each point gets 1/3 of the total area
    weights = np.array([area / 3.0, area / 3.0, area / 3.0])

    return points, weights



#6.1
def poissons_lhs(u, v, grad_u, grad_v, x, dx):
    return dot(grad_u, grad_v)*dx
    

#6.2
def poisson_rhs(v, grad_v, x, dx):
    a=0.25
    return exp(-((x[0]- 0.5)**2 +(x[1]-0.5)**2)/(2*a**2))*v*dx


#6.3
def compute_element_matrix_2d(lhs, X):
    Ak = np.zeros((3, 3))
    #kvadratur punkter och vikter
    points, weights = generate_quadrature_2d(X)

    for k in range(len(weights)):
        x = points[:, k]
        dx = weights[k]
        
        for i in range(3):
            v, grad_v = evaluate_basis_2d(i, x, X)
            for j in range(3):
                u, grad_u = evaluate_basis_2d(j, x, X)
                Ak[i, j] = Ak[i, j] + lhs(u, v, grad_u, grad_v, x, dx)
    return Ak


#6.4
def compute_element_vector_2d(rhs, X):
    Bk = np.zeros(3)

    points, weights = generate_quadrature_2d(X)
    for k in range(len(weights)):
        x = points[:, k]
        dx = weights[k]
        for i in range(3):
            v, grad_v = evaluate_basis_2d(i, x, X)
            Bk[i] += rhs(v, grad_v, x, dx)
    
    return Bk

#6.5
from numpy import array
from scipy.sparse import lil_matrix


def assemble_matrix_2d(lhs, V, K):

    N = len(V)
    A = lil_matrix((N, N))

    
    for K_ in K:
        X = array([V[K_[0]], V[K_[1]], V[K_[2]]]).T

        A_K = compute_element_matrix_2d(lhs, X)
        # Add element matrix to stiffness matrix
        for i in range(3):
            I = K_[i]
            for j in range(3):
                J = K_[j]
                A[I, J] += A_K[i, j]
    print(A)
    return A


#6.6
def assemble_vector_2d(rhs, V, K):
    N = len(V)
    b = zeros(N)


    for K_ in K:
        X = array([V[K_[0]], V[K_[1]], V[K_[2]]]).T

   
        b_K = compute_element_vector_2d(rhs, X)

        
        for i in range(3):
            I = K_[i]
            b[I] += b_K[i]

    return b

#6.7


def apply_dirichlet_bc_2d(A, b, V):

    # Create zero row
    N = shape(V)[0]
    zero = zeros(N)

    # Iterate over vertices
    for i, v in enumerate(V):

        # Check if we are on the boundary
        eps = 1e-6
        if v[0] < eps or v[0] > 1 - eps or v[1] < eps or v[1] > 1 - eps:

            # Zero row
            A[i, :] = zero

            # Insert 1 on the diagonal
            A[i, i] = 1

            # Set Dirichlet value in vector
            b[i] = 0


#6.8
from scipy.sparse.linalg import spsolve



def solve_poisson_2d(V, K):

    # Assemble stiffness matrix
    A = assemble_matrix_2d(poissons_lhs, V, K)

    # Assemble load vector
    b = assemble_vector_2d(poisson_rhs, V, K)

    # Apply boundary conditions
    apply_dirichlet_bc_2d(A, b, V)

    # Compute solution
    U = spsolve(A.tocsr(), b)

    return U



            


