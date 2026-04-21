import numpy as np
from scipy.linalg import eig
from mtm026 import assem, Ke_balk, Ke_sigma_balk, extract_block

# Parametrar
P0 = 1.0; L1 = 5.0; L2 = 5.0
a = 0.05; t = 0.002; E = 200e9

# Tröghetsmomenter
Iv = (a**4 - (a - 2*t)**4) / 12
Ih = (a * (2*a)**3 - (a - 4*t) * (2*a - 4*t)**3) / 12
EIv = E * Iv
EIh = E * Ih

K = np.zeros((14, 14))
K_sig = np.zeros((14, 14))

Ke_v = Ke_balk(EIv, L2)
Ke_h = Ke_balk(EIh, L1)
Ke_sig_v = Ke_sigma_balk(P0, L2)

assem(K, Ke_v, dofs=[1, 2, 5, 6])
assem(K_sig, Ke_sig_v, dofs=[1, 2, 5, 6])

assem(K, Ke_v, dofs=[5, 6, 10, 11])
assem(K_sig, Ke_sig_v, dofs=[5, 6, 10, 11])

assem(K, Ke_v, dofs=[3, 4, 5, 8])

assem(K, Ke_v, dofs=[5, 8, 10, 13])

assem(K, Ke_h, dofs=[7, 6, 9, 8])

assem(K, Ke_h, dofs=[12, 11, 14, 13])

free_dofs = [3, 4, 5, 6, 5, 8, 10, 11, 10, 13]

k_red = extract_block(K, free_dofs, free_dofs)
k_red_sig = extract_block(K_sig, free_dofs, free_dofs)

alpha , vec = eig(k_red, k_red_sig)

print(vec[3])
print(np.min(alpha)*P0)