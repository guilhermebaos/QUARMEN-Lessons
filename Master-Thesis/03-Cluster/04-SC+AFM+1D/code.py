# Usual stuff
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# System
import sys

# Typing
from typing import Callable

# Path
from pathlib import Path


# Read inputs
g_input = float(sys.argv[1])
u_input = float(sys.argv[2])


# Directory where the current script is located
script_path = Path(__file__).resolve()
script_name = script_path.stem
script_dir = script_path.parent



# GENERAL CODE -----

def bissect(func: Callable, a: float, b: float, eps: float = 1e-6, maxI: int = 200, mult: bool = False) -> float:
    """
        Solves the equation `func(x) = 0` using the bissection method in [a, b].

        The function stops when the error is smaller then `eps` or the number of iterations exceeds `maxI`.

        The error is additive by default, but can be made multiplicative by setting `mult` = True

        ### Return
        x: The value of x that solves the equation.
    """

    # Evaluate on the edges
    fa, fb = func(a), func(b)

    # Check if there is a zero on the edges
    if abs(fa) <= 1e-16:
        return a
    elif abs(fb) <= 1e-16:
        return b

    # Main Loop
    i = 0
    error = (b-a)/2

    # For multiplicative error make sure to do at least one iteration
    if mult:
        error = eps + 1

    while abs(error) > eps:

        # Find the midpoint and evaluate the function there
        c = (a+b)/2
        fc = func(c)

        # Check if the zero is exactly on the midpoint
        if abs(fc) <= 1e-16:
            return c

        # Sign change happens in [a, c]
        elif fa * fc < 0:
            b = c
            fb = fc
        
        # Sign change happens in [c, b]
        elif fc * fb < 0:
            a = c
            fa = fc
        
        # No sign change
        else:
            raise ValueError(f"O sinal de f é o mesmo em a = {a} e em b = {b} (f(a) = {fa}, f(b) = {fb}) - iteração {i}).")

        # Compute the error
        error = (b - a)/2

        # If the error is multiplicative, compare the error with the midpoint value
        if mult:
            error /= (a+b)/2

        # Maximum number of iterations
        i += 1
        if i > maxI:
            raise ValueError(f"Não foi possível encontrar a raíz da função com a precisão desejada em menos de maxI = {maxI} iterações!")
        
    return (a+b)/2


def eps(kk: np.ndarray) -> np.ndarray:
    """
        Compute the dispersion relation for a given set of vectors `kk`.
    """
    return -2 * np.sum(np.cos(kk), axis = 0)


def fock(decimal: int, modes: int = 8) -> list[int]:
    return [(decimal >> n) % 2 for n in range(modes)][::-1]


# Sample code to convert from state to full fock space
# The four least significant bits (LSB) corresponds to k and the four MSBs corresponds to k+pi
state = 20
n_pk_ppi_u_i, n_pk_ppi_d_i, n_mk_ppi_u_i, n_mk_ppi_d_i, n_pk_npi_u_i, n_pk_npi_d_i, n_mk_npi_u_i, n_mk_npi_d_i = fock(state, 8)



# CONSTANTS
MODES = 8
DIM = 2**MODES


# SUBSPACES
class Subspace:

    # Get properties and states
    # Technically we could even forget about the quantum numbers
    def __init__(self, parity: int, spin: int) -> None:
        # Quantum numbers of the subspace
        self.parity = parity
        self.spin = spin

        # States and dimension
        self.states = []
        self.dimension = 0
        self.Nk = 0

        # Lookup table
        self.lookup = dict()
    
    # Print
    def __str__(self) -> str:
        return f"Space with P = {self.parity} and S = {self.spin} contains {len(self.states)} states."

    
    # STATES
    
    # Add a state to the list
    def add_state(self, st: int) -> None:
        self.states += [st]

    # Build a lookup table {state: position} on the list
    def build_lookup(self) -> None:
        self.dimension = len(self.states)

        for index, st in enumerate(self.states):
            self.lookup[st] = index


    # OPERATORS

    # Create the operators
    def build_operators(self) -> tuple[np.ndarray, np.ndarray]:

        # Make sure we have a lookup table
        if len(self.lookup) == 0: self.build_lookup()

        # Operators (t for total, meaning both options)
        self.nk_tk_tpi_t = np.zeros((self.dimension, self.dimension), dtype=np.float64)

        self.nk_tk_tpi_d = np.zeros((self.dimension, self.dimension), dtype=np.float64)
        self.nk_tk_tpi_u = np.zeros((self.dimension, self.dimension), dtype=np.float64)

        self.nk_tk_npi_t = np.zeros((self.dimension, self.dimension), dtype=np.float64)
        self.nk_tk_ppi_t = np.zeros((self.dimension, self.dimension), dtype=np.float64)
        
        self.bk_tk_tpi = np.zeros((self.dimension, self.dimension), dtype=np.float64)

        self.afmk = np.zeros((self.dimension, self.dimension), dtype=np.float64)

        for st in self.states:
            # Position of the state
            ket = self.lookup[st]

            # Binary decomposition
            st_list = fock(st)
            n_pk_ppi_u, n_pk_ppi_d, n_mk_ppi_u, n_mk_ppi_d, n_pk_npi_u, n_pk_npi_d, n_mk_npi_u, n_mk_npi_d = st_list

            # Number of particles in the state
            self.nk_tk_tpi_t[ket, ket] = np.sum(st_list)

            # Number of particles with spin up and spin down
            self.nk_tk_tpi_u[ket, ket] = n_pk_ppi_u + n_mk_ppi_u + n_pk_npi_u + n_mk_npi_u
            self.nk_tk_tpi_d[ket, ket] = n_pk_ppi_d + n_mk_ppi_d + n_pk_npi_d + n_mk_npi_d

            # Number of particles at ±k and ±(k+pi)
            self.nk_tk_npi_t[ket, ket] = n_pk_npi_u + n_pk_npi_d + n_mk_npi_u + n_mk_npi_d
            self.nk_tk_ppi_t[ket, ket] = n_pk_ppi_u + n_pk_ppi_d + n_mk_ppi_u + n_mk_ppi_d


            # Applying bk to this state (st is the ket, we find the bra)
            # For the first and third lines, the phase is n_pk_npi_u-1 + n_pk_npi_d + n_mk_npi_u bu the first term is 0 because it is 1-1
            if n_mk_npi_d * n_pk_npi_u:
                bra_pk_npi = self.lookup[st - 8 - 1]
                self.bk_tk_tpi[bra_pk_npi, ket] += (-1)**(n_pk_npi_d + n_mk_npi_u)
            
            if n_pk_npi_d * n_mk_npi_u:
                bra_mk_npi = self.lookup[st - 4 - 2]
                self.bk_tk_tpi[bra_mk_npi, ket] += (-1)**(n_pk_npi_d)
            
            if n_pk_ppi_u * n_mk_ppi_d:
                bra_pk_ppi = self.lookup[st - 128 - 16]
                self.bk_tk_tpi[bra_pk_ppi, ket] += (-1)**(n_pk_ppi_d + n_mk_ppi_u)
            
            if n_pk_ppi_d * n_mk_ppi_u:
                bra_mk_ppi = self.lookup[st - 64 - 32]
                self.bk_tk_tpi[bra_mk_ppi, ket] += (-1)**(n_pk_ppi_d)
            
            # Staggered magnetization measurement
            if (1 - n_pk_npi_u) * n_pk_ppi_u:
                bra_pk_u = self.lookup[st + 8 - 128]
                self.afmk[bra_pk_u, ket] += (-1)**(n_pk_ppi_d + n_mk_ppi_u + n_mk_ppi_d)
                
            if (1 - n_pk_npi_d) * n_pk_ppi_d:
                bra_pk_d = self.lookup[st + 4 - 64]
                self.afmk[bra_pk_d, ket] += -(-1)**(n_mk_ppi_u + n_mk_ppi_d + n_pk_npi_u)

            if (1 - n_mk_npi_u) * n_mk_ppi_u:
                bra_mk_u = self.lookup[st + 2 - 32]
                self.afmk[bra_mk_u, ket] += (-1)**(n_mk_ppi_d + n_pk_npi_u + n_pk_npi_d)
                
            if (1 - n_mk_npi_d) * n_mk_ppi_d:
                bra_mk_d = self.lookup[st + 1 - 16]
                self.afmk[bra_mk_d, ket] += -(-1)**(n_pk_npi_u + n_pk_npi_d + n_mk_npi_u)
            
        # Symmetrize AFM (it is an observable)
        self.afmk += self.afmk.T
            
        return self.nk_tk_tpi_t, self.bk_tk_tpi
            
        

    def build_hamiltonian(self, Nk: int, ee_npi: np.ndarray, ee_ppi: np.ndarray, mu: float, U: float, Delta: complex, J: float, m: float) -> np.ndarray:

        # Number of k-points
        self.Nk = Nk

        # Dispersion relation
        xi_npi = ee_npi - mu
        xi_ppi = ee_ppi - mu

        # Hamiltonian
        self.ham = np.zeros((Nk, self.dimension, self.dimension), dtype=np.float64)

        # Diagonal
        for st in self.states:
            # Position of the state
            ket = self.lookup[st]

            # Binary decomposition
            st_list = fock(st)
            n_pk_ppi_u, n_pk_ppi_d, n_mk_ppi_u, n_mk_ppi_d, n_pk_npi_u, n_pk_npi_d, n_mk_npi_u, n_mk_npi_d = st_list

            # Kinetic Energy
            self.ham[:, ket, ket] += xi_npi * np.sum(n_pk_npi_u + n_pk_npi_d + n_mk_npi_u + n_mk_npi_d)
            self.ham[:, ket, ket] += xi_ppi * np.sum(n_pk_ppi_u + n_pk_ppi_d + n_mk_ppi_u + n_mk_ppi_d)

            # HK Interaction
            self.ham[:, ket, ket] += U * (n_pk_ppi_u * n_pk_ppi_d + n_mk_ppi_u * n_mk_ppi_d + n_pk_npi_u * n_pk_npi_d + n_mk_npi_u * n_mk_npi_d)
            
            # AFM Interaction (nn method)
            # self.ham[:, ket, ket] += -J * (n_pk_npi_u * n_pk_ppi_d + n_mk_npi_u * n_mk_ppi_d + n_pk_npi_d * n_pk_ppi_u + n_mk_npi_d * n_mk_ppi_u)
            # self.ham[:, ket, ket] += -J * ((1 - n_pk_npi_u) * (1 - n_pk_ppi_d) + (1 - n_mk_npi_u) * (1 - n_mk_ppi_d) + (1 - n_pk_npi_d) * (1 - n_pk_ppi_u) + (1 - n_mk_npi_d) * (1 - n_mk_ppi_u))

            # self.ham[:, ket, ket] += -J * (n_pk_npi_u * n_pk_ppi_u + n_mk_npi_u * n_mk_ppi_u + n_pk_npi_d * n_pk_ppi_d + n_mk_npi_d * n_mk_ppi_d)

            # self.ham[:, ket, ket] += -J * ((n_pk_npi_u - n_pk_npi_d) * (n_pk_ppi_u - n_pk_ppi_d) + (n_mk_npi_u - n_mk_npi_d) * (n_mk_ppi_u - n_mk_ppi_d))

            # self.ham[:, ket, ket] += -J * (n_pk_npi_u * n_pk_ppi_u + n_pk_npi_d * n_pk_ppi_d - n_pk_npi_u * n_pk_ppi_d - n_pk_npi_d * n_pk_ppi_u)
            # self.ham[:, ket, ket] += -J * (n_mk_npi_u * n_mk_ppi_u + n_mk_npi_d * n_mk_ppi_d - n_mk_npi_u * n_mk_ppi_d - n_mk_npi_d * n_mk_ppi_u)

            # FM Interaction (for Testing)
            # self.ham[:, ket, ket] += -J * 10 * (n_pk_ppi_u + n_mk_ppi_u + n_pk_npi_u + n_mk_npi_u)
        

        # Off-diagonal

        # SC Interaction
        self.ham += np.conjugate(Delta) * self.bk_tk_tpi
        self.ham += Delta * np.transpose(self.bk_tk_tpi)

        # AFM Interaction (Mean-Field Method)
        # TODO: Figure out why this extra factor of 2 is needed to align with the theory!
        self.ham += J * 2 * m * self.afmk

        # with np.printoptions(precision=3, suppress=True, threshold=100_000, linewidth=200):
        #     print(f"The space is (P, S) = ({self.parity}, {self.spin})")

        return self.ham



def generate_fspace() -> dict[tuple[int, int], Subspace]:
    """
        Generate the Fock space for our problem.
    """

    # Create all spaces (types so we can have hints)
    fspace = dict[tuple[int, int], Subspace]()
    for parity in [0, 1]:
        for spin in range(-MODES//2, MODES//2+1):
            # If we have an even number of particles then spin must be even
            if abs(spin) % 2 == parity:
                fspace[(parity, spin)] = Subspace(parity, spin)


    # STATES
    for st in range(DIM):

        # Compute parity
        parity = np.sum([(st >> n) % 2 for n in range(MODES)]) % 2

        # Compute spin
        spin = np.sum([(st >> (2*n + 1)) % 2 for n in range(MODES//2)]) - np.sum([(st >> (2*n)) % 2 for n in range(MODES//2)])

        # Add state to correct subspace
        fspace[(parity, spin)].add_state(st)


    # Build lookup tables
    # Check size of spaces, and that they add to the whole space
    maxi = 0
    total = 0
    for key in fspace:
        fspace[key].build_lookup()
        maxi = np.max([maxi, len(fspace[key].states)])
        total += len(fspace[key].states)

        # print(key, len(fspace[key].states))

    return fspace

    # print(f"Is the total number of states {DIM}? {"Yes" if total == DIM else "No"}")
    # print(f"The biggest space has {maxi} states.")



def thermal_average(fspace: dict[tuple[int, int], Subspace], eigen: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]], opcode: str, T: float, E0: np.ndarray, mult_kk: np.ndarray = np.array([]), make_positive: bool = False) -> float:

    # Number of k-points
    Nk = fspace[(0, 0)].Nk

    # Partition function
    Z = np.zeros(Nk)

    # Boltzmann weights
    bb = dict()


    for key in fspace:
        # Get eigenvalues and eigenvecs
        vals, _ = eigen[key]

        # Minimum energy, to scale Boltzmann factors
        valsE0 = vals - E0[:, None]

        # Boltzmann exponential
        bexp = np.exp(-valsE0/T)
        bb[key] = bexp

        # Compute partition function (sum over all dim states for each k)
        Z += np.sum(bexp, axis = 1)

    # Compute Boltzmann factor
    for key in eigen:

        # Divide by partition function
        bb[key] /= Z[:, None]
    
    # Compute average via a trace
    res = 0
    for key in fspace:
        
        # Get eigenvalues and eigenvecs
        vals, vecs = eigen[key]
        
        # Choose operator
        op = getattr(fspace[key], opcode)
        
        # Choose Boltzman factor
        bfac = bb[key]

        # Operator times the kets
        temp = op @ vecs
        
        # Bra times (op |ket>) via expectation[k, i] = sum_j conj(vecs[k, j, i]) * temp[k, j, i]
        expectation = np.einsum("kji,kji->ki", vecs.conj(), temp)

        if make_positive:
            expectation = np.abs(expectation)

        # Sum over i with Boltzmann factor
        res_kk = np.sum(bfac * expectation, axis=1)

        if len(mult_kk) == 0:
            # Sum over k
            res += np.sum(res_kk)
        
        else:
            # Sum over k with the factor mult for each k-point
            res += np.sum(res_kk * mult_kk)

    return res


def solve(L: int, d: int, mu: float, U: float, g: float, J: float, T: float, delta_start: float = 1e-2, delta_eps: float = 1e-4, m_start: float = 1e-2, m_eps: float = 1e-4, alpha: float = 0.8):
    """
        Solve the mean-field Hamiltonian for `L` k-points along each `d` directions for temperature `T` and chemical potential `mu`.

        Self-consistently compute Delta starting at `delta_start` with an error of `delta_eps`.
    """

    # Possible adjustments
    # Take out the [:-1] from kk
    # Change the kx > 0 into kx >= 0


    # OBTAINING K-POINTS

    # Determine values of k in our lattice
    kk = np.linspace(-np.pi, np.pi, L + 1)[:-1]

    # Make a mesh of k values
    grid_kk = np.meshgrid(*[kk]*d, indexing='ij')

    # Convert to a single array with shape (d, L ** d)
    mesh_kk = np.stack(grid_kk).reshape(d, len(kk) ** d)

    # Select only 0 <= kx <= np.pi/2
    mesh_kk = mesh_kk[:, (0 < mesh_kk[0]) *  (mesh_kk[0] <= np.pi/2)]

    # Compute corresponding energies
    ee_npi = eps(mesh_kk)
    ee_ppi = eps(mesh_kk + np.pi)

    # Number of k-points
    Nk = len(mesh_kk[0])


    # SET UP FOCK SPACE

    # Generate the space
    fspace = generate_fspace()

    # Prepare the subspaces
    for key in fspace:
        fspace[key].build_operators()


    # SELF-CONSISTENT CALCULATION OF DELTA
    if g == 0:
        delta_start = 0

    if J == 0:
        m_start = 0

    delta_error = delta_eps + 1
    m_error = m_eps + 1
    while delta_error > delta_eps or m_error > m_eps:
        
        # Solve the Hamiltonian
        E0 = np.zeros(Nk) + 1e6
        eigen = dict()
        for key in fspace:

            # Generate Hamiltonian
            ham = fspace[key].build_hamiltonian(Nk, ee_npi, ee_ppi, mu, U, delta_start, J, m_start)

            # Diagonalize Hamiltonians
            vals, vecs = np.linalg.eigh(ham)

            # Save the results
            eigen[key] = (vals, vecs)
        
            # Compute minimum energy for each k-point
            E0 = np.minimum(E0, np.min(vals, axis=1))
            
        # Compute the component of delta for this subspace
        delta_new = thermal_average(fspace, eigen, "bk_tk_tpi", T, E0)
        delta_new *= -g / (4 * Nk)

        # Compute error
        delta_error = np.abs(delta_start - delta_new)

        # Prepare next loop
        delta_start = alpha * delta_new + (1 - alpha) * delta_start


        # Compute the component of m for this subspace
        m_new = thermal_average(fspace, eigen, "afmk", T, E0)
        m_new *= -1 / (2 * 4 * Nk)

        # Compute error
        m_error = np.abs(m_start - m_new)

        # Prepare next loop
        m_start = alpha * m_new + (1 - alpha) * m_start
    

    # SOLVE FOR CONVERGED ORDER PARAMETERS
    E0 = np.zeros(Nk) + 1e6
    eigen = dict()
    for key in fspace:

        # Generate Hamiltonian
        ham = fspace[key].build_hamiltonian(Nk, ee_npi, ee_ppi, mu, U, delta_start, J, m_start)
        
        # Diagonalize Hamiltonians
        vals, vecs = np.linalg.eigh(ham)

        # Save the results
        eigen[key] = (vals, vecs)
    
        # Compute minimum energy for each k-point
        E0 = np.minimum(E0, np.min(vals, axis=1))


    # Compute expectation values
    delta = thermal_average(fspace, eigen, "bk_tk_tpi", T, E0)
    n = thermal_average(fspace, eigen, "nk_tk_tpi_t", T, E0)

    afm = thermal_average(fspace, eigen, "afmk", T, E0)
    fm = thermal_average(fspace, eigen, "nk_tk_tpi_u", T, E0) - thermal_average(fspace, eigen, "nk_tk_tpi_d", T, E0)

    # Kinetic energy has different contributions from npi and ppi
    Kxx = 0
    Kxx += thermal_average(fspace, eigen, "nk_tk_npi_t", T, E0, mult_kk=-ee_npi)
    Kxx += thermal_average(fspace, eigen, "nk_tk_ppi_t", T, E0, mult_kk=-ee_ppi)

    # Normalize
    delta *= -g / (4 * Nk)
    n *= 1 / (4 * Nk)
    afm *= -1 / (2 * 4 * Nk)
    fm *= 1 / (4 * Nk)

    Kxx *= np.pi / (4 * Nk)

    return np.abs(delta), n, afm, fm, Kxx



# TASK-SPECIFIC CODE ----- TODO


# PARAMETERS (CHANGE THE VALUES)


# CONSTANTS
d = 1
W = 4 * d

# PARAMETERS
L = 800
g = g_input * W
U = u_input * W
T = 1e-9 * W

# CALCULATION PARAMETERS
# We use higher errors while computing the chemical potential mu

# Delta
delta_start = 1e-2
delta_eps_mu = 1e-4
delta_eps = 1e-6

# Magnetization
m_start = 1e-2
m_eps_mu = 1e-4
m_eps = 1e-5


# J Sweep
J_min = 0*W
J_max = 0.15*W
J_ste = 13

jj = np.linspace(J_min, J_max, J_ste)


# X sweep
x_min = 0
x_max = 0.3
x_ste = 13

xx = np.linspace(x_min, x_max, x_ste)
mu_eps = 4e-4


# Calculations

# Make the mesh
xx_mesh, jj_mesh = np.meshgrid(xx, jj)

dd_mesh = np.zeros_like(jj_mesh)
mm_mesh = np.zeros_like(jj_mesh)
for ij, J in enumerate(jj):

    print(f"Step {ij + 1}/{J_ste}: computing for J / W = {J / W:.2f}")

    succ = 0
    for ix, x in enumerate(xx):

        # Target filling
        n_target = 1 - x

        # Bissect for mu
        mu_min = -W/2 - U - g - J
        mu_max = U/2 + W

        mu = bissect(lambda mu_test: solve(L, d, mu_test, U, g, J, T, delta_start=delta_start, delta_eps=delta_eps_mu, m_start=m_start, m_eps=m_eps_mu)[1] - n_target, mu_min, mu_max, eps = mu_eps)
        
        # Solve for converged mu
        delta, n, fm, _, _ = solve(L, d, mu, U, g, J, T, delta_start=delta_start, delta_eps=delta_eps, m_start=m_start, m_eps=m_eps)

        # Save
        dd_mesh[ij, ix] = delta
        mm_mesh[ij, ix] = fm

        # Show progress
        print(f"\rProgress: {ix + 1}/{x_ste}", end="")

    # Reset print
    print("")


# OUTPUT
script_out = script_dir / "outputs"
script_out.mkdir(exist_ok=True)

np.savez(script_out / f"{d}={L}={T}={g}={U}", xx_mesh=xx_mesh, jj_mesh=jj_mesh, dd_mesh=dd_mesh, mm_mesh=mm_mesh)
