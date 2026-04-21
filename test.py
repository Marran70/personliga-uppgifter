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

# Elementmatriser
Ke_v     = np.array(Ke_balk(EIv,  L2)).astype(float)
Ke_h     = np.array(Ke_balk(EIh,  L1)).astype(float)
Ke_sig_v = np.array(Ke_sigma_balk(P0, L2)).astype(float)

# Bundna DOF (platshållare):
#   8  = u_x vänster botten = 0  (fast inspänd)
#   9  = u_y vänster botten = 0  (fast inspänd)
#   10 = θ   vänster botten = 0  (fast inspänd)
#   11 = u_y höger botten   = 0  (rullstöd, endast u_y låst)
# Fria DOF:
#   1 = θ    vänster mitt
#   2 = θ    höger mitt
#   3 = u_x  mitt-nivå (delad)
#   4 = θ    vänster topp
#   5 = θ    höger topp
#   6 = u_x  topp-nivå (delad)
#   7 = θ    höger botten  ← FRI (rullstöd tillåter rotation)

ndof = 11
K     = np.zeros((ndof, ndof))
K_sig = np.zeros((ndof, ndof))

# Elem 1: vänster vertikal botten→mitt
# [v_bot, θ_bot, v_mitt, θ_mitt] = [u_x=0, θ=0, u_x_mitt, θ_mittV]
#                                 = [8,     10,   3,        1      ]
assem(K,     Ke_v,     dofs=[8, 10, 3, 1])
assem(K_sig, Ke_sig_v, dofs=[8, 10, 3, 1])

# Elem 2: vänster vertikal mitt→topp
assem(K,     Ke_v,     dofs=[3, 1, 6, 4])
assem(K_sig, Ke_sig_v, dofs=[3, 1, 6, 4])

# Elem 3: höger vertikal botten→mitt
# [v_bot, θ_bot, v_mitt, θ_mitt] = [u_x=0, θ_fri, u_x_mitt, θ_mittH]
#                                 = [8,     7,      3,         2      ]
# OBS: u_x höger botten = 0 pga axialdef=0 → samma bundna DOF 8 som vänster
assem(K, Ke_v, dofs=[8, 7, 3, 2])

# Elem 4: höger vertikal mitt→topp
assem(K, Ke_v, dofs=[3, 2, 6, 5])

# Elem 5: horisontell balk mitt-nivå
# v (vertikal utböjning) ≈ 0 i båda ändar → platshållare DOF 11
assem(K, Ke_h, dofs=[11, 1, 11, 2])

# Elem 6: horisontell balk topp-nivå
assem(K, Ke_h, dofs=[11, 4, 11, 5])

# Extrahera 7×7 reducerad matris
free_dofs = [1, 2, 3, 4, 5, 6, 7]
K_red     = extract_block(K,     free_dofs, free_dofs)
K_sig_red = extract_block(K_sig, free_dofs, free_dofs)

# Lös egenvärdesproblemet
vals, vecs = eig(K_red, K_sig_red)
real_vals  = np.real(vals)
alpha_kr   = np.min(real_vals[real_vals > 0])
P_kr       = alpha_kr * P0

print(f"Lägsta lastmultiplikator  α_kr = {alpha_kr:.2f}")
print(f"Kritisk knäcklast         P_kr = {P_kr/1000:.4f} kN")

# Knäckningsmod
mode_idx = np.where(np.isclose(real_vals, alpha_kr))[0][0]
mode     = np.real(vecs[:, mode_idx])
dof_names = ["θ_mittV", "θ_mittH", "u_x_mitt", "θ_toppV", "θ_toppH", "u_x_topp", "θ_högerBot"]
print("\nEgenvektorn:")
for name, val in zip(dof_names, mode):
    print(f"  {name:12s} = {val:.4f}")