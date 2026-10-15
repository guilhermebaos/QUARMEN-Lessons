using SymPy

# 1. Define the symbolic parameters
U, t_AA, t_BB = symbols("U t_AA t_BB", real=true)
t_AB = symbols("t_AB")
t_BA = conj(t_AB)

println("--- Analytical First Eigenvalue ---")
println("lambda_0 = ", t_AA + t_BB)
println("v_0 = [0, 1/sqrt(2), 0, -1/sqrt(2)]\n")

# 2. Construct the reduced 3x3 matrix in the { e1, e+, e3 } basis
# Using SymPy's S(2) to ensure exact square roots
sq2 = sqrt(sympy.S(2))

M_reduced = [
    2*t_AA + U      -sq2*t_AB        0;
    -sq2*t_BA        t_AA + t_BB     -sq2*t_AB;
    0               -sq2*t_BA        2*t_BB + U
]

# Convert to native SymPy Matrix object
M3 = sympy.Matrix(M_reduced)

println("--- Reduced 3x3 Matrix ---")
display(M3)

# 3. Find the Eigenvalues of the 3x3 block
println("\n--- Eigenvalues of the 3x3 block ---")
evals = M3.eigenvals()
display(evals)

# 4. Diagonalize the 3x3 Matrix
println("\n--- Diagonalizing 3x3 Matrix ---")
P_3, D_3 = M3.diagonalize()

println("\nTransformation Matrix P_3 (Eigenvectors in reduced basis):")
display(P_3)

println("\nDiagonalized Matrix D_3 (Eigenvalues on diagonal):")
display(D_3)