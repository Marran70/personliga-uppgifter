import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from hjälpfunktioner import generate_mesh_2d, assemble_matrix_2d
#nu ska vi upprepda exakt samma algoritm egentligen men nu finns en ny del i ekvationen eftersom V = 10^3(x1-1/2)^2
#detta gör att vi får ännu en lhs funktion att beräkna och denna kallas potential. 
#denna ser nästan ut som massmatrisen som multiplicerar funktionsvärdena men här också med V alltså 
#integralen över V(x1, x2) * v * u*dx i svag form vilket gör att vi får den nya potential funktionen nedan. 

#fulltändiga ekvationen på svag form är alltså: 
#integral av 1/2 grad_U * grad_v*dx + integral av V(x1, x2)*u*v*dx = lambda*integral av u*v*dx

#Steg 1: Definiera lhs-funktioner
#Kinetisk del, samma som tidigare
def kinetic_lhs(u, v, grad_u, grad_v, x, dx):
    return 0.5 * np.dot(grad_u, grad_v) * dx

#Potentialdel V(x1, x2)*u*v, V utvärderas i kvadraturpunkten x
def potential_lhs(u, v, grad_u, grad_v, x, dx):
    V = 1e3 * (x[0] - 0.5)**2   #V(x1, x2) = 10^3*(x1 - 0.5)^2
    return V * u * v * dx        #ingen gradient, som massmatrisen

#Massmatris, samma som tidigare
def mass_lhs(u, v, grad_u, grad_v, x, dx):
    return u * v * dx
# Steg 2: Assemblera A = A_kinetisk + A_potential

def bygg_och_los(n):
    V_noder, K = generate_mesh_2d(n, n)

    # Assemblera de två bidragen till A separat och addera
    A_kin = assemble_matrix_2d(kinetic_lhs,   V_noder, K)
    A_pot = assemble_matrix_2d(potential_lhs, V_noder, K)
    A = A_kin + A_pot

    M = assemble_matrix_2d(mass_lhs, V_noder, K)
    #Hitta inre noder (samma som tidigare)
    eps = 1e-10
    interior = []
    for i in range(len((V_noder))):
        v = V_noder[i]
        if v[0] > eps and v[0] < 1-eps and v[1] > eps and v[1] < 1-eps:
            interior.append(i)

    #Extrahera sub-matriser ( samma som tidigare)
    A_int = A[interior, :][:, interior].tocsr()
    M_int = M[interior, :][:, interior].tocsr()

    return V_noder, K, interior, A_int, M_int

# Steg 3: Plotta 4 lägsta egenfunktionerna (32x32-nät)

n = 32
V_noder, K, interior, A_int, M_int = bygg_och_los(n)
num_eigs = 4
eigenvalues, eigenvectors = eigsh(A_int, k=num_eigs, M=M_int, which='SM')

# Sortera minst först
idx = np.argsort(eigenvalues)
eigenvalues  = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Mappa tillbaka till alla noder (samma som i (a))
N = len(V_noder)
U_full = np.zeros((N, num_eigs))
for k in range(num_eigs):
    U_full[interior, k] = eigenvectors[:, k]

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle('Fyra lägsta egenfunktioner med V = 10^3(x_1-1/2)^2', fontsize=13)

#samma for loop som i (a)
for k in range(num_eigs):
    ax = axes[k // 2, k % 2]
    tcf = ax.tricontourf(V_noder[:, 0], V_noder[:, 1], K, U_full[:, k], levels=20, cmap='RdBu_r')
    fig.colorbar(tcf, ax=ax)
    ax.set_title(f'lambda_{k+1} = {eigenvalues[k]:.4f}')
    ax.set_xlabel('x_1')
    ax.set_ylabel('x_2')

plt.tight_layout()
#plt.savefig('stjärn_egenfunktioner.png', dpi=150)
plt.show()

# Steg 4:Konvergensstudie – jämför mot finaste nätet

# Ingen analytisk lösning finns, så vi använder 128x128 som referens
print("Beräknar referensvärde på 128x128-nät:")
_, _, interior_ref, A_ref, M_ref = bygg_och_los(128)
lambda_ref, _ = eigsh(A_ref, k=1, M=M_ref, which='SM') #beräkna egenvärden på samma sätt som tidigare, k=1 pga behöver bara det minsta
lambda_ref = lambda_ref[0]  #ta detta som referensvärde
print(f"Referensvärde lambda_ref = {lambda_ref:.8f}")

# Loopa över grövre nät
n_values = [8, 16, 32, 64]
h_values = []
errors   = []

for n in n_values:  #loopa över nät på samma sätt som i uppgift (b) och ska ska vi sedan jämföra dessa grövre med referensen. 
    print(f"Beräknar för {n}x{n}-nät:")
    _, _, interior_n, A_n, M_n = bygg_och_los(n)
    lam, _ = eigsh(A_n, k=1, M=M_n, which='SM') #återigen bara minsta som behövs
    error = abs(lam[0] - lambda_ref)  #beräknar error med hjälp av vårat referensvärde
    h_values.append(1/n)
    errors.append(error) #spara de olika errorsen för de olika näten. 
    print(f" h = {1/n:.4f},  lambda_h = {lam[0]:.6f},  fel = {error:.2e}")

# Konvergensordning (samma sätt som i uppgift (b))
log_h = np.log(h_values)
log_e = np.log(errors)
p = np.polyfit(log_h, log_e, 1)[0]
print(f"\nKonvergensordning p = {p:.2f}")

# Plotta
fig, ax = plt.subplots(figsize=(7, 5))
ax.loglog(h_values, errors, 'bo-', label='FEM-fel lambda_h - lambda_ref')
h_ref = np.array(h_values)
ref_line = errors[0] * (h_ref / h_values[0])**p
ax.loglog(h_ref, ref_line, 'r--', label=f'Referens h^{p:.1f}')
ax.set_xlabel('h (nätstorlek)', fontsize=12)
ax.set_ylabel('Fel', fontsize=12)
ax.set_title('Konvergens med V = 10^3(x_1-1/2)^2', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
#plt.savefig('stjärn_konvergens.png', dpi=150)
plt.show()