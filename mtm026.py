import sympy as sp
import numpy as np
from scipy.linalg import eig, eigh
from IPython.display import Math, display
from numbers import Real
from typing import List

"""
Hjälpfunktioner vid problemlösning i kursen Hållfasthetslära - MTM026
Jim Brouzoulis, 25-03-2025

Korrigerad av Amanda Hansen 27-03-2026
"""


def extract_block(K, rows:List[int], cols:List[int]=[0])-> np.ndarray:
    """
    Extrahera en blockmatris från en matris `K` genom att ange rader `rows` och kolumner `cols` som ska plockas ut.
    Frihetsgraders numrering börjar på 1.

    Indata
    ----------
    K : np.ndarray
        Den globala styvhetsmatrisen
    rows:
        lista av vilka rader som ska axtraheras
    cols:
        lista av vilka kollumner som ska axtraheras

    Utdata
    -------
    np.ndarray
        Den extraherade delen av styvhetsmatrisen som svarar mot dessa rader och kollumner

    Exempel
    -------
    >>> K = np.array([[10, 2, 3],
    ...               [2, 20, 5],
    ...               [3, 5, 30]])
    >>> rows = [1, 3]
    >>> cols = [1, 3]
    >>> extract_block(K, rows, cols)
    array([[10,  3],
           [ 3, 30]])
    """
    nrows = len(rows)
    ncols = len(cols)
    if issubclass(type(K), sp.Matrix): # om K är en symbolisk matris returnera en symbolisk matris
        K_red = sp.zeros(nrows, ncols)
    else: # antag numpy array
        K_red = np.zeros((nrows, ncols), dtype=K.dtype)
    for i, row in enumerate(rows):
        for j, col in enumerate(cols):
            K_red[i, j] = K[row-1, col-1]
    return K_red

def assem(K:np.ndarray, Ke:np.ndarray, dofs:List[int]):
    """Assemplering av matriser

    Indata:
        K (np.ndarray): global K-matris (styvhetsmatris)
        Ke (np.ndarray): elementets styvhetsmatris
        dofs (List[int]): elementets frihetsgrader, frihetsgraders numrering börjar på 1.

    Utdata:
        K (np.ndarray): Uppdaterad global styvhetsmatris

    Exempel:

        assem(K=K, Ke=Ke1, dofs=[1,2,3,4])

    Notering:
        Denna funktion går att använda till vektorer också, exempelvis kraftvektorn.

    """
    # Kontroll av input
    nrows, ncols = K.shape
    # kontrollera att frihetsgraderna finns i matrisen
    if np.max(dofs) > nrows:
        raise AssertionError(f"Du försöker assemblera till frihetsgrader som inte finns")
    # assert np.max(dofs) <= nrows

    # Kontrollera att ingen frihetsgrad har index 0
    if np.min(dofs) <= 0:
            raise AssertionError(f"Alla frihetsgrader måste ha nummer större än 0. dofs = {dofs}")
    # assert np.min(dofs) > 0 # kontrollera att ingen frihetsgrad har index 0

    if nrows == ncols: # kvadratisk matris
    # assert nrows == ncols # kontrollera att det är en kvadratisk matris
    # assert np.max(dofs) <= nrows - 1 # kontrollera att frihetsgraderna finns i matrisen
        for row, dof_i in enumerate(dofs):
            for col, dof_j in enumerate(dofs):
                K[dof_i-1, dof_j-1] += Ke[row, col]
    elif ncols == 1: # vektor
        for row, dof_i in enumerate(dofs):
            K[dof_i-1, 0] += Ke[row]


    return K


def Ke_balk(EI:Real, L:Real):
    """
    Styvhetsmatris för ett balkelement med fyra frihetsgrader [a₁ a₂ a₃ a₄], där a₁ och a₃ är utböjning samt a₂ och a₄ är rotationer.

    Indata:
        * `EI` - Böjstyvhet, produkten mellan E och I. Antas vara konstant över elementet.
        * `L` - Balkelementets längd

    Utdata:
        * `Ke` - styvhetsmatris, storlek [4×4]

    Exempel:

        Ke = Ke_balk(EI=1, L=2)

    """
    return  EI / L**3 * sp.Matrix([
        [12, 6*L, -12, 6*L],
        [6*L, 4*L**2, -6*L, 2*L**2],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*L**2, -6*L, 4*L**2]
    ])

def Ke_sigma_balk(P:Real, L:Real):
    """
    Spännings-styvhetsmatris för ett balkelement med fyra frihetsgrader [a₁ a₂ a₃ a₄], där a₁ och a₃ är utböjning samt a₂ och a₄ är rotationer.

    Indata:
        * `P` - Tryckkraft (N(x) = - P)
        * `L` - Balkelementets längd

    Utdata:
        * `Ke_σ ` - styvhetsmatris, storlek [4×4]

    Exempel:

        Ke = Ke_sigma_balk(P=1, L=2)

    """
    return  P / (30*L) * sp.Matrix([
        [36, 3*L, -36, 3*L],
        [3*L, 4*L**2, -3*L, -L**2],
        [-36, -3*L, 36, -3*L],
        [3*L, -L**2, -3*L, 4*L**2]
    ])

def fe_balk(q:Real, L:Real):
    """
    Ekvivalent nodlastvektor för balkelement med längd `L` belastad med en utbredd last med konstant intensitet `q`
    """
    return sp.Matrix([
        [q*L/2],
        [q*L**2/12],
        [q*L/2],
        [-q*L**2/12]
        ])

def fe_balk_linjär_last(q:Real, L:Real):
    """
    Ekvivalent nodlastvektor för balkelement med längd `L` belastad med en
    utbredd last med linjärt minskande intensitet från `q` till 0, från nod 1-2
    """
    return sp.Matrix([
        [7*q*L/20],
        [q*L**2/20],
        [3*q*L/20],
        [-q*L**2/30]
        ])

def fe_balk_linjär(q1:Real, q2:Real, L:Real):
    """
    Ekvivalent nodlastvektor för balkelement med längd `L` belastad med en
    utbredd last med linjärt varierande intensitet från q₁ i nod 1 till q₂ i nod 2
    """
    return np.array([
        [(7*q1 + 3*q2)*L/20],
        [L**2*(q1/20 + q2/30)],
        [(3*q1 + 7*q2)*L/20],
        [-L**2*(q1/30 + q2/20)]
        ])

def Ke_fjäder(k:Real):
    """
    Styvhetsmatris för ett fjäderelement med styvhet k med två frihetsgrader [a₁ a₂], där a₁ och a₂ kan vara axiella förskjutningar
    i de två noderna eller rotationer i de tvä ändarna för en rotationsfjäder.

    Indata:
        * `k` - Fjäderstyvhet, antas vara konstant över elementet.
        * `k = E*A/L` - för en stång

    Utdata:
        * `Ke` - styvhetsmatris, storlek [2×2]

    Exempel:

        Ke = Ke_fjäder(k=1)

    """
    return  k * sp.Matrix([
        [1, -1],
        [-1, 1],
    ])

def displayvar(name:str, var, accuracy:int=None):
    """
    Skriv ut en variabel `var` med variabelnamnet `name` på formen: name = var
    accuracy - avrunda till decimalform med givet antal värdesiffror. Default är att skriva ut de exakta uttrycken.
    Exempel:

        displayvar("P", 1)

    """
    if isinstance(var, np.ndarray):
        var = sp.Matrix(var)
    if accuracy is None:
        display(Math(f'{name} = {sp.latex(var)}'))
    else:
        display(Math(f'{name} \\approx {sp.latex(sp.sympify(var).evalf(accuracy))}'))


def display_eqnsys(K, a, f, accuracy:int=None):
    """
    Skriv ut ett ekvationssystem på formen: K a = f
    accuracy - avrunda till decimalform med givet antal värdesiffror. Default är att skriva ut de exakta uttrycken.
    """
    if accuracy is None:
        display( Math( f"{sp.latex(sp.Matrix(K))} {sp.latex(a)} = {sp.latex(f)}" ) )
    else:
        display( Math( f"{sp.latex(sp.Matrix(K).evalf(accuracy))} {sp.latex(a)} \\approx {sp.latex(f.evalf(accuracy))}" ) )
