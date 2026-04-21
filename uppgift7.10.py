from mtm026 import *

P0 = L = EI = 1

# med låg fjäderkonstant leder detta tilla tt vi får 
#under = E2
#över = E3

K = np.zeros((5, 5))
K_sig = np.zeros((5, 5))
Ke = Ke_balk(EI=EI, L=L)
Ke_sig = Ke_sigma_balk(P = P0, L=L)
assem(K, Ke, dofs=[1, 2, 3, 4])
assem(K_sig, Ke_sig, dofs=[1, 2, 3, 4])

k_fjäder=Ke_fjäder(k = 1)
assem(K, k_fjäder, dofs=[4, 5])

free_dofs = [2, 4]
K_red = extract_block(K, free_dofs, free_dofs)
K_sig_red = extract_block(K_sig, free_dofs, free_dofs)
alpha, a = eig(K_red, K_sig_red)

Alphrak = np.min(alpha)
print(Alphrak*P0)
