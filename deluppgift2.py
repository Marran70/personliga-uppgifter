import numpy as np
from scipy.linalg import eig
from mtm026 import assem, Ke_balk, Ke_sigma_balk, extract_block

# Parametrar
P0 = 1.0; L1 = 5.0; L2 = 5.0; a = 0.05; t = 0.002; E = 200e9

# Beräkna styvheter (EI)
Iv = (a**4 - (a - 2*t)**4) / 12
Ih = (a * (2*a)**3 - (a - 4*t) * (2*a - 4*t)**3) / 12
EIv = E*Iv 
EIh = E*Ih

# Definiera matriser (t.ex. 10x10 för att få plats med alla index)
K = np.zeros((10, 10))
K_sig = np.zeros((10, 10))

# Element-matriser
Ke_v = Ke_balk(EIv, L2)
Ke_h = Ke_balk(EIh, L1)
Ke_sig_v = Ke_sigma_balk(P0, L2)

# Fria frihetsgrader:
# 1: u_x mitt, horisontell förskjutning av mittvåningen
#  2: theta_mitt_V, rotation i vänster nod mitt
#  3: theta_mitt_H, rotation i höger nod mitt
#  4: u_x topp,   horisontell förskutning av toppvåning
# 5: theta_topp_V,  rotation i vänster nod topp
#  6: theta_topp_H  rotation i höger nod topp
# Låsta index (platshållare): 7, 8, 9, 10

# Assemblera
# Vertikala vänstra (bär last)
# Element 1: vänster vertikal botten→mitt
# Ke_balk DOF: [v_botV, θ_botV, v_mittV, θ_mittV]
# Bundna: v_botV=0 (DOF 7), θ_botV=0 (DOF 8)
# Fria:   v_mittV = u_x_mitt (DOF 1), θ_mittV (DOF 2)
assem(K,     Ke_v,     dofs=[7, 8, 1, 2])
assem(K_sig, Ke_sig_v, dofs=[7, 8, 1, 2])

# Element 2: vänster vertikal mitt→topp
assem(K,     Ke_v,     dofs=[1, 2, 4, 5])
assem(K_sig, Ke_sig_v, dofs=[1, 2, 4, 5])

# Element 3: höger vertikal botten→mitt
assem(K, Ke_v, dofs=[9, 10, 1, 3])

# Element 4: höger vertikal mitt→topp
assem(K, Ke_v, dofs=[1, 3, 4, 6])

# Element 5: horisontell balk mittnivå  ← ÄNDRAT
assem(K, Ke_h, dofs=[1, 2, 1, 3])

# Element 6: horisontell balk toppnivå  ← ÄNDRAT
assem(K, Ke_h, dofs=[4, 5, 4, 6])

# Extrahera de fria delarna
free_dofs = [1, 2, 3, 4, 5, 6]
K_red = extract_block(K, free_dofs, free_dofs)
K_sig_red = extract_block(K_sig, free_dofs, free_dofs)

# Lös egenvärdesproblemet
vals, vecs = eig(K_red, K_sig_red)
real_vals = np.real(vals)
alpha_kr  = np.min(real_vals[real_vals > 0])
P_kr = alpha_kr * P0
print(f"Kritisk last P_kr: {P_kr/1000:.2f} kN")