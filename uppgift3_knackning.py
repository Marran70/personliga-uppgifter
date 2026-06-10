from mtm026 import *


# Parametrar
E   = 200e9
P0  = 1.0
L   = 5.0   # L1 = L2 = 5 m
a   = 0.05
t   = 0.002

# Tröghetsmomenter
# Vertikala balkar: kvadratiskt ihåligt tvärsnitt a x a, vägg t
I_vert  = (a**4 - (a - 2*t)**4) / 12

# Horisontella balkar: rektangulärt ihåligt tvärsnitt a x 2a, vägg 2t
I_horis = (a * (2*a)**3 - (a - 4*t) * (2*a - 4*t)**3) / 12

EI_v = E * I_vert
EI_h = E * I_horis

# -------------------------------------------------------
# Reducerad modell – 6 frihetsgrader efter att axialdef. försummats:
#   q1 = u_C = u_D  (horis. förflyttning mittplan, länkas via horis. balk)
#   q2 = θ_C        (rotation vänster mitt)
#   q3 = θ_D        (rotation höger mitt)
#   q4 = u_E = u_F  (horis. förflyttning toppplan)
#   q5 = θ_E        (rotation vänster topp)
#   q6 = θ_F        (rotation höger topp)
#
# Inlåsta: θ_A = θ_B = u_A = u_B = 0 (inspänd botten) → dof = 0 i assem()
# Vertikala förflyttningar = 0 (axialdef. försummas i vertikala balkar)
# -------------------------------------------------------

K     = np.zeros((6, 6))
K_sig = np.zeros((6, 6))

# Hjälpfunktion: Ke_balk/Ke_sigma_balk returnerar sympy-matriser → konvertera
def np_(M): return np.array(M).astype(float)

# Element 1: Vänster undre pelare (A→C), tryckkraft P0
# [v1,θ1, v2,θ2] → [u_A,θ_A, u_C,θ_C], inspänning sätter u_A=θ_A=0 → dofs=[1,2]
assem(K,     np_(Ke_balk(EI=EI_v, L=L)),        dofs=[1, 2])
assem(K_sig, np_(Ke_sigma_balk(P=P0, L=L)),      dofs=[1, 2])

# Element 2: Vänster övre pelare (C→E), tryckkraft P0
assem(K,     np_(Ke_balk(EI=EI_v, L=L)),         dofs=[1, 2, 4, 5])
assem(K_sig, np_(Ke_sigma_balk(P=P0, L=L)),       dofs=[1, 2, 4, 5])

# Element 3: Höger undre pelare (B→D), tryckkraft = 0 (lasten verkar bara vänster)
# inspänning sätter u_B=θ_B=0 → dofs=[1,3]
assem(K, np_(Ke_balk(EI=EI_v, L=L)), dofs=[1, 3])

# Element 4: Höger övre pelare (D→F), tryckkraft = 0
assem(K, np_(Ke_balk(EI=EI_v, L=L)), dofs=[1, 3, 4, 6])

# Element 5: Undre horisontell balk (C→D)
# Axialdef. försummas i pelarna → w_C = w_D = 0 → bidrar bara med rotationsstyvhet
# [v1,θ1, v2,θ2] → [0, θ_C, 0, θ_D] → dofs=[2,3]
assem(K, np_(Ke_balk(EI=EI_h, L=L)), dofs=[2, 3])

# Element 6: Övre horisontell balk (E→F), samma resonemang → dofs=[5,6]
assem(K, np_(Ke_balk(EI=EI_h, L=L)), dofs=[5, 6])

# Lös generaliserat egenvärdesproblem: (K_red - α * K_sig_red) * a = 0
free_dofs = [1, 2, 3, 4, 5, 6]
K_red     = extract_block(K,     free_dofs, free_dofs)
K_sig_red = extract_block(K_sig, free_dofs, free_dofs)

alpha, a_vec = eig(K_red, K_sig_red)

# Välj reella, positiva egenvärden och sortera
mask   = (np.abs(alpha.imag) < 1e-6 * (np.abs(alpha.real) + 1e-30)) & (alpha.real > 0)
alphas = np.sort(alpha[mask].real)

print(f"Lastmultiplikatorer α = {alphas[:4]}")  # välj den lägsta
print(f"Kritisk last P_kr = {alphas[0] * P0:.2f} N = {alphas[0] * P0 / 1000:.4f} kN")

# P_kr ≈ 18 590 N = 18.59 kN  → uppfyller kravet P_kr > 15 kN (del 3.2)
