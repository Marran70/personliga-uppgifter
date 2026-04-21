#undre E2
#övre E3

from mtm026 import *
P0 = L = EI = 1

K = np.zeros((6, 6))
K_sig = np.zeros((6, 6))

Dofs = [
    [1, 2, 3, 4], 
    [3, 4, 5, 6]
]

for i in range(2):
    dofs = Dofs[i]
    K1 = Ke_balk(EI, L=L)
    
    assem(K, K1, dofs)


K1_sig=Ke_sigma_balk(P=P0, L=L) #denna inte i loop pga P bara verkar i mitten av balken(se bild)
assem(K_sig, K1_sig, dofs)

free_dofs = [2, 4, 6]

K_red = extract_block(K, free_dofs, free_dofs)
K_red_sig= extract_block(K_sig, free_dofs, free_dofs)

alpha, a = eig(K_red, K_red_sig)

alphak = np.min(alpha)
print(alphak*P0)