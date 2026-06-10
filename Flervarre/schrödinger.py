import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from hjälpfunktioner import generate_mesh_2d, assemble_matrix_2d, compute_element_matrix_2d
#Steg 1: Skapa nätet

#generate_mesh_2d(32,32) skapar (32+1)*(32+1) = 1089 noder
L = 1
nx, ny = 32, 32
V, K = generate_mesh_2d(nx, ny)

#Steg 2: Definiera lhs-funktioner
#Styvhetsmatris A: (1/2) * grad_u · grad_v
def schrodinger_lhs(u, v, grad_u, grad_v, x, dx):
    return 0.5 * np.dot(grad_u, grad_v) * dx

# Massmatris M: u * v  använder bara funktionsvärden, ingen gradient
#vi måste ha denna då lösningen av Au=lam*u antar att integralen av u*v = 1 vilket inte är sant
#så vi måste istället lösa Au = lam*M*u vilket är det vi gör sen. 
def mass_lhs(u, v, grad_u, grad_v, x, dx):
    return u * v * dx
# Steg 3: Assemblera matriserna
# assemble_matrix_2d loopar över alla element K och
# anropar compute_element_matrix_2d för varje triangel
A = assemble_matrix_2d(schrodinger_lhs, V, K)
M = assemble_matrix_2d(mass_lhs, V, K)

# Steg 4: Hitta inre noder (ej på randen)
# Randnoder har x=0, x=1, y=0 eller y=1
# Vi tar bara med inre noder i egenvärdesproblemet, det är det denna for loop hittar
eps = 1e-10
interior = []
for i in range(len((V))):
        v = V[i]
        if v[0] > eps and v[0] < 1-eps and v[1] > eps and v[1] < 1-eps:
            interior.append(i)
# Steg 5: Extrahera sub-matriser för inre noder

# Vi plockar ut raderna och kolumnerna som svarar mot inre noder
# Detta är ekvivalent med att sätta u=0 på randen
A_int = A[interior, :][:, interior].tocsr() #styvhetsmatrsien för endast inre noder
M_int = M[interior, :][:, interior].tocsr() #Massmatrisen men bara för inre noder
#.tocsr gör den till en sparse-format om vi kan lösa egennvärdespoblemet med

# Steg 6: Lös generaliserat egenvärdes­problem
# eigsh hittar de k minsta egenvärden för A*u = lambda*M*u
# which='SM' = Smallest Magnitude
num_eigs = 4
eigenvalues, eigenvectors = eigsh(A_int, k=num_eigs, M=M_int, which='SM')

# Sortera efter storleksordning (minst först)
idx = np.argsort(eigenvalues)
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Mappa tillbaka till alla noder (randnoder = 0)
N = len(V)
U_full = np.zeros((N, num_eigs))
for k in range(num_eigs):
    U_full[interior, k] = eigenvectors[:, k]

# Steg 7: Plotta egenfunktionerna med tricontourf

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle('De fyra lägsta egenfunktionerna (L=1)', fontsize=14)

for k in range(num_eigs):
    ax = axes[k // 2, k % 2]
    # tricontourf(x-kord, y-kord, trianglar, värden)
    tcf = ax.tricontourf(V[:, 0], V[:, 1], K, U_full[:, k], levels=20, cmap='RdBu_r')
    fig.colorbar(tcf, ax=ax)
    ax.set_title(f'lambda{k+1} = {eigenvalues[k]:.4f}') #sätt titel för varje graf
    ax.set_xlabel('x1') #vilket axel det visar
    ax.set_ylabel('x2')

plt.tight_layout()
plt.savefig('egenfunktioner.png', dpi=150)
plt.show()

# Jämför med analytiska egenvärden

print("Numeriska egenvärden:")  #vilka värden vi får när vi räknar ut med python
for k in range(num_eigs):
    print(f"  lambda_{k+1} = {eigenvalues[k]:.6f}")

print("\nAnalytiska egenvärden lambda_mn = (n^2+m^2)*pi^2/2:")

analytiska = []  #vilka värden vi fick på papper i tidigare uppgiften. 
for n in range(1, 4):
     for m in range(1, 4):
          lam = (n**2 + m**2) * np.pi**2 / 2
          analytiska.append(lam)
analytiska = sorted(analytiska)

for k in range(len(analytiska[:6])):
    lam = analytiska[k]
    print(f"   lambda_{k+1} = {lam:.6f}")