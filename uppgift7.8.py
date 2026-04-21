from mtm026 import *

#undre gräns: euler 1 pga om fjädern minskar så blir strukturen vekare som om det inte finns någon fjäder
#viker då enligt euler 1 som är fri i toppen. 
#övre gräns: Detta blir en euler 3 pga det till slut inte kommer kunna röra sig i sidled utan bara upp och ner, samt är fast inspänd i botten

p0= E = I = L = 1
EI = E*L


P_E1 = (np.pi**2)*EI/(4*L^2)
P_E3 = (2.05*(np.pi**2))*EI/(L**2)

K = np.zeros((5, 5))
K_sig = np.zeros((5, 5))
K1 = Ke_balk(EI=EI, L=L)
K1_sig = Ke_sigma_balk(P=p0, L=L)
assem(K, K1, dofs=[1, 2, 3, 4])
assem(K_sig, K1_sig, dofs=[1, 2, 3, 4])

K_fjäder = Ke_fjäder(k=1)
assem(K, K_fjäder, dofs=[3, 5])

free_dofs = [3, 4]
K_red = extract_block(K, free_dofs, free_dofs)
k_sigma_red = extract_block(K_sig, free_dofs, free_dofs)
alpha , a = eigh(K_red, k_sigma_red)

alphakr = np.min(alpha)
print(alphakr)

print(alphakr*p0)


