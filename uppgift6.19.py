from mtm026 import *

EI, P, L, ra, rb, rc = sp.symbols("EI, P, L, Ra, Rb, Rc")
q0 = P/L
Dofs = [
    [1, 2, 3, 4],
    [3, 4, 5, 6],
    [5, 6, 7, 8]
]

a1, a2, a3, a4, a5, a6, a7, a8 = sp.symbols("a1:9")
f1, f2, f3,f4, f5, f6, f7, f8 = sp.symbols("f1:9")

a = sp.Matrix([0, a2, 0, a4, a5, a6, 0, a8])

fb = sp.Matrix([ra, 0, rb, 0, -P, 0, rc, 0])

K = sp.zeros(8, 8)
fl = sp.zeros(8, 1)
fle = fe_balk(q=-q0, L=L/3)
for el in range(3):
    dofs = Dofs[el]
    
    Ke = Ke_balk(EI=EI, L=L/3)
    assem(fl, fle, dofs)
    assem(K, Ke, dofs)

obekanta = [a2, a4, a5, a6, a8, ra, rb, rc]
ekv_sys = sp.Eq(K*a, fb + fl)
sol = sp.solve(ekv_sys, obekanta)

display(ekv_sys)
display(a.subs(sol))
display(fl.subs(sol))