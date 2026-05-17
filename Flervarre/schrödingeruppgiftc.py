import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from hjälpfunktioner import generate_mesh_2d, assemble_matrix_2d
#hela detta programmet fungerar egentligen typ exakt samma som det första bara att vi loopad igenom de olika näten
#samt att vi sedan beräknas konvergensordningen samt plottar den på en log log graf 
#och jämför med en refenslinje vilket är den linjen som den vi beräkna bör följa. 

#Steg 1: Definiera lhs-funktioner (samma som i (b))
def schrodinger_lhs(u, v, grad_u, grad_v, x, dx):
    return 0.5 * np.dot(grad_u, grad_v) * dx

def mass_lhs(u, v, grad_u, grad_v, x, dx):
    return u * v * dx

#Steg 2: Analytiskt egenvärde att jämföra mot
lambda_exact = np.pi**2  #lambda_11 = (1^2+1^2)pi^2/2 = pi^2

#Steg 3: Loopa över nätsstorlekar

#n_values är antalet element i varje riktning
#h = 1/n är steglängden
n_values = [8, 16, 32, 64, 128]
h_values = []
errors = []

for n in n_values:
    print(f"Beräknar för {n}x{n}-nät:")
    #Skapa nätet för aktuella n_values
    h = 1 / n
    V, K = generate_mesh_2d(n, n)
    #Assemblera matriserna för aktuellt nät
    A = assemble_matrix_2d(schrodinger_lhs, V, K)
    M = assemble_matrix_2d(mass_lhs, V, K)

    #Hitta inre noder (samma som i (b))
    eps = 1e-10
    interior = []
    for i in range(len((V))):
        v = V[i]
        if v[0] > eps and v[0] < 1-eps and v[1] > eps and v[1] < 1-eps:
            interior.append(i)

    #Extrahera sub-matriser för inre noder ( samma som i (b))
    A_int = A[interior, :][:, interior].tocsr()
    M_int = M[interior, :][:, interior].tocsr()

    #Beräkna minsta egenvärdet
    #k=1 eftersom vi bara behöver det minsta, det är det som vi ska jämföra med sen. 
    eigenvalues, _ = eigsh(A_int, k=1, M=M_int, which='SM')
    lambda_h = eigenvalues[0]

    #Beräkna felet
    error = abs(lambda_h - lambda_exact)
    
    #Spara h och felet
    h_values.append(h)
    errors.append(error)

    print(f" h = {h:.4f},  lambda_h = {lambda_h:.6f},  fel = {error:.2e}")

#Steg 4: Beräkna konvergensordningen
#I en log-log plot ska felet vara en rät linje: log(errors) = p*log(h) + C
#Vi skattar lutningen p med numpy.polyfit
log_h = np.log(h_values)
log_e = np.log(errors)

# polyfit(x, y, 1) passar en rät linje till datan
# koefficienter[0] är lutningen = konvergensordningen p
koefficienter = np.polyfit(log_h, log_e, 1)
p = koefficienter[0]
print(f"\nKonvergensordning p = {p:.2f}")

# Steg 5: Plotta log-log

fig, ax = plt.subplots(figsize=(7, 5))
# Plotta felet
ax.loglog(h_values, errors, 'bo-', label='FEM-fel lambda_h - lambda_11')

# Plotta en referenslinje med lutning p för jämförelse
# Vi skalar den så att den passar vår data
h_ref = np.array(h_values)
ref_line = errors[0] * (h_ref / h_values[0])**p
ax.loglog(h_ref, ref_line, 'r--', label=f'Referens h^{p:.1f}') #referenslinje

ax.set_xlabel('h: (nätstorlek)', fontsize=12)
ax.set_ylabel('Fel: lambda_h - lambda_11', fontsize=12)
ax.set_title('Konvergens för minsta egenvärdet', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
#plt.savefig('konvergens.png', dpi=150)
plt.show()