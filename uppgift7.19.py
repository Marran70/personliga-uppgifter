#undre E3
#övre E4

from mtm026 import *
P0 = L = EI = 1

K = np.zeros((8, 8))
K_sig = np.zeros((8, 8))

Dofs = [
    [1, 2, 3, 4], 
    [3, 4, 5, 6],
    [5, 6, 7, 8]
]
el_length = [L/2, L/2, L]
for i in range(3):
    dofs = Dofs[i]
    K1 = Ke_balk(EI, L=el_length[i])
    
    assem(K, K1, dofs)


K1_sig=Ke_sigma_balk(P=P0, L=L/2) 
assem(K_sig, K1_sig, dofs=[1, 2, 3, 4])
assem(K_sig, K1_sig, dofs=[3, 4, 5, 6])

free_dofs = [3, 4, 6, 8]

K_red = extract_block(K, free_dofs, free_dofs)
K_red_sig= extract_block(K_sig, free_dofs, free_dofs)

alpha, a = eig(K_red, K_red_sig)

print(alpha)
print(np.min(alpha)*P0)