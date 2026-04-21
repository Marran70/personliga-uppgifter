#gränserna är E1 på båda men med olika längd samt olika EI

from mtm026 import *

P0 = EI = L = 1

Dofs = [
    [1,2, 3, 4],
    [3, 4, 5, 6]
]
K = np.zeros((6, 6))
K_sig = np.zeros((6, 6,))
Shuno = [2*EI, EI]
for i in range(len(Dofs)):
    dofs = Dofs[i]
    ke = Ke_balk(EI=Shuno[i], L=L)
    ke_sig=Ke_sigma_balk(P=P0, L=L)
    assem(K, ke, dofs)
    assem(K_sig, ke_sig, dofs)

free_dofs = [3, 4, 5, 6]
k_red = extract_block(K, free_dofs, free_dofs)
k_red_sig = extract_block(K_sig, free_dofs, free_dofs)

alpha, a = eigh(k_red, k_red_sig)

print(alpha)

print(a)

print(np.min(alpha)*P0)
