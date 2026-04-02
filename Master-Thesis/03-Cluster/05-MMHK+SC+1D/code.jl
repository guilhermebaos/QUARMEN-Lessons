# Imports
import Base


# ----- HELPER FUNCTIONS -----
"""
    bissect(func, a, b, eps, maxI, mult)

Find a root of the function `func` within the interval `[a, b]` using the bisection method.

### Arguments
- `func::Function`: The continuous function to evaluate.
- `a::Real`, `b::Real`: The lower and upper bounds of the search interval.
- `eps::Real`: The error threshold used as stopping criterion.
- `maxI::Int`: The maximum number of iterations allowed before throwing an error.
- `mult::Bool`: If `true`, uses a multiplicative (relative) error; if `false`, uses an additive (absolute) error.

### Returns
- `Real`: An approximate solution to f(x) = 0 with error `eps`.

### Errors
- `ArgumentError`: Thrown if no sign change is detected within the sub-interval or if `maxI` is exceeded.
"""
function bissect(func::Function, a::Real, b::Real, eps::Real, maxI::Int, mult::Bool)::Real

    # Evaluate on the edges
    fa, fb = func(a), func(b)

    # Check if there is a zero on the edges
    if abs(fa) <= 1e-16
        return a
    elseif abs(fb) <= 1e-16
        return b
    end

    # Main Loop
    i = 0

    # Starting error (make sure to do at least one iteration)
    error = eps + 1

    while abs(error) > eps

        # Find the midpoint and evaluate the function there
        c = (a + b) / 2
        fc = func(c)

        # Check if the zero is exactly on the midpoint
        if abs(fc) <= 1e-16
            return c

        # Sign change happens in [a, c]
        elseif fa * fc < 0
            b = c
            fb = fc
        
        # Sign change happens in [c, b]
        elseif fc * fb < 0
            a = c
            fa = fc
        
        # No sign change
        else
            throw(ArgumentError("The sign of f is the same at a = $a and at b = $b because f(a) = $fa, f(b) = $fb - Iteration $i."))
        end

        # Compute the error
        error = (b - a)/2

        # If the error is multiplicative, compare the error with the midpoint value
        # We must have a safeguard against division by zero, which may occur if the root is at zero
        if mult
            mid = (a + b) / 2
            error /= abs(mid) < 1e-14 ? 1.0 : mid
        end

        # Maximum number of iterations
        i += 1
        if i > maxI
            throw(ArgumentError("The root of function f cannot be found to the desired precision in less then maxI = $maxI iterations"))
        end
    end

    return (a + b) / 2
end


"""
    disp(kk::AbstractArray)

Compute the dispersion relation for a given set of vectors `kk` which are *in 1D*.
"""
function disp(kk::AbstractArray)
    return -2 * cos.(kk)
end


"""
    fock(decimal::Int, modes::Int=8)

Convert a decimal integer into a binary vector corresponding to its Fock state representation.

Each entry represents whether the mode is occupied (1) or empty (0),
meaning the bitstring has length `modes`.

# Arguments
- `decimal::Int`: The integer value to convert.
- `modes::Int`: The number of modes (bits) to include.

# Returns
- `Vector{Int}`: A list of bits representing the Fock state.
"""
function fock(decimal::Int, modes::Int)::Vector{Int}
    return digits(decimal, base=2, pad=modes)[end:-1:1]
end

function fock!(decimal::Int, modes::Int, output::Vector{Int})
    output .= digits(decimal, base=2, pad=modes)[end:-1:1]
end


# TODO: Replace fock with a non-allocating version using bit shifts

#=
FOCK SPACE

Our Fock space will contain 2^(4n) states, due to +k mixing with -k and each k-point having 2n states
where 2 comes from the 2-element fermionic Fock basis per state and 4n is the number of states (2 k-points * 2 spins * n orbitals).

We organize it as follows, for any number of mixed momenta n, in order of decreasing priority:
1. First the +k states then the -k states.
2. First the lower-indexed orbitals.
3. First the up-spin then the down-spin.

For n = 1 we have (where n(k, α, σ) is the occupation number):
|n(+k, 1, ↑), n(+k, 1, ↓), n(-k, 1, ↑), n(-k, 1, ↓)>

And for n = 2:
|n(+k, 1, ↑), n(+k, 1, ↓), n(+k, 2, ↑), n(+k, 2, ↓), n(-k, 1, ↑), n(-k, 1, ↓), n(-k, 2, ↑), n(-k, 2, ↓)>
=#

"""

# Arguments
- `k::Int`: Use k = 0 for +k and k = 1 for -k.
- `alpha::Int`: Use 1, ..., n to select the orbital.
- `sigma::Int`: Use sigma = 0 for ↑ and sigma = 1 for ↓

# Returns
- `Int`: The position of the desired mode on the Fock state according to our convention
"""
function npos(k::Int, alpha::Int, sigma::Int, nMMHK::Int)
    return 1 + k * (2 * nMMHK) + (alpha - 1) * 2 + sigma
end


# ----- SUBSPACES -----
struct propSubspace
    # Subspace quantum numbers
    spin::Int64

    # States in the subspace
    states::Vector{Int64}

    # Dimension of the space
    dimension::Int64

    # Lookup table for states
    lookup::Dict{Int64, Int64}
end


struct opsSubspace
    # Subspace properties
    prop::propSubspace

    # Operator list
    # Parameterization allows Int matrices for nk operators and Complex for the Hamiltonian
    ops::Dict{String, Matrix{Int64}}
end


struct hamSubspace
    # Subspace properties
    prop::propSubspace

    # Operator list
    # Parameterization allows Int matrices for nk operators and Complex for the Hamiltonian
    ops::Dict{String, Matrix{Int64}}

    # Hamiltonian
    ham::Matrix{ComplexF32}
end


# Print the subspace
function Base.show(io::IO, sub::propSubspace)
    print(io, "Space with S = $(sub.spin) contains $(length(sub.states)) states: $(sub.states).")
end


# Build all subspaces
# We work with a temporary dict to build the spaces, so that they can be immutable
function buildPropSubs(modes::Int)
    # Dimension of the Fock space
    dim = 2^modes

    # Possible spins
    fspace_ss = collect(-div(modes, 2):div(modes, 2))
    max_ss = div(modes, 2)

    # Go through all states and add them to the correct spin subspace
    fspace_states = [Int64[] for _ in 0:modes]
    fspace_lookup = [Dict{Int64, Int64}() for _ in 0:modes]
    for state in 0:dim-1
        # Compute spin
        spin = sum(n * (2 * (pos % 2) - 1) for (pos, n) in enumerate(fock(state, modes)))

        # Compute position on the list
        spin_index = spin + max_ss + 1

        # Store in the proper list
        push!(fspace_states[spin_index], state)

        # Build lookup table
        # Keys are states as integers and values are their position on the list
        fspace_lookup[spin_index][state] = length(fspace_states[spin_index])
    end


    # Check that the states add to the total number of states
    total = sum(length(x) for x in fspace_states)
    if total != dim
        throw(ErrorException("Subspaces are incorrect!"))
    end

    # Build the subspaces
    fspace = [propSubspace(spin, states, length(states), lookup) for (spin, states, lookup) in zip(fspace_ss, fspace_states, fspace_lookup)]

    return fspace
end






# ----- OPERATORS -----
function buildOpsSubs(sub::propSubspace, nMMHK::Int, modes::Int)
    # Create the operators
    ops = Dict{String, Matrix{Int64}}()

    # Empty operators
    ops["nk_tk_ts"] = zeros(Int64, sub.dimension, sub.dimension)

    ops["nk_tk_up"] = zeros(Int64, sub.dimension, sub.dimension)
    ops["nk_tk_dw"] = zeros(Int64, sub.dimension, sub.dimension)
    
    ops["nk_pk_ts"] = zeros(Int64, sub.dimension, sub.dimension)
    ops["nk_mk_ts"] = zeros(Int64, sub.dimension, sub.dimension)

    ops["tk"] = zeros(Int64, sub.dimension, sub.dimension)
    
    ops["bk_tk"] = zeros(Int64, sub.dimension, sub.dimension)
    ops["fmk"] = zeros(Int64, sub.dimension, sub.dimension)


    # Fill in the operators
    st_list = fock(0, modes)
    for st in sub.states
        # Position of the state
        ket = sub.lookup[st]

        # Binary decomposition, by mutating st_list
        fock!(st, modes, st_list)

        # Number of particles
        ops["nk_tk_ts"][ket, ket] = sum(st_list)

        # Number of particles with a given spin (odds are spin up, evens are spin down)
        ops["nk_tk_up"][ket, ket] = sum(n * (pos % 2) for (pos, n) in enumerate(st_list))
        ops["nk_tk_dw"][ket, ket] = sum(n * ((pos + 1) % 2) for (pos, n) in enumerate(st_list))

        # Number of particles at +k and -k
        ops["nk_pk_ts"][ket, ket] = sum(n * (pos <= modes / 2) for (pos, n) in enumerate(st_list))
        ops["nk_mk_ts"][ket, ket] = sum(n * (pos >  modes / 2) for (pos, n) in enumerate(st_list))


        # Applying bk to this state for each alpha (st is the ket, we find the bra)
        # The operator is given by bk = c(-p, α, ↓)c(p, α, ↑)
        for alpha in 1:nMMHK
            # Applying bk for +k
            # if nk_mk_ds * nk_pk_ds for this orbital
            pos_pk_us = npos(0, alpha, 0, nMMHK)
            pos_mk_ds = npos(1, alpha, 1, nMMHK)
            if (st_list[pos_pk_us] * st_list[pos_mk_ds]) != 0
                # The bra is the state without those two particles
                bra_pk = sub.lookup[st - 2^(modes - pos_pk_us) - 2^(modes - pos_mk_ds)]

                # The phase is given by the states between the destroyed particles
                ops["bk_tk"][bra_pk, ket] = (-1)^sum(st_list[(pos_pk_us + 1):(pos_mk_ds - 1)]; init=0)
            end

            # Applying bk for -k
            # if nk_pk_ds * nk_mk_ds for this orbital
            pos_pk_ds = npos(0, alpha, 1, nMMHK)
            pos_mk_us = npos(1, alpha, 0, nMMHK)
            if (st_list[pos_pk_ds] * st_list[pos_mk_us]) != 0
                # The bra is the state without those two particles (the -1 is because Julia is 1-indexed)
                bra_mk = sub.lookup[st - 2^(modes - pos_pk_ds) - 2^(modes - pos_mk_us)]

                # The phase is given by the states between the destroyed particles
                # The +1 is from the fact that now we destroy the -k first, which has to commute with the filled +k state
                ops["bk_tk"][bra_mk, ket] = (-1)^(sum(st_list[(pos_pk_ds + 1):(pos_mk_us - 1)]; init=0) + 1)
            end
        end
    end

    # Magnetization
    # TODO: Make magnetization per orbital
    ops["fmk"] .= ops["nk_tk_up"] - ops["nk_tk_dw"]

    # Return the subspace
    return opsSubspace(sub, ops)
end


"""
    Create the Hamiltonian, without any of the self-consistent parameters!
"""
function buildHamSubs(sub::opsSubspace, nMMHK::Int, modes::Int, k::Real, U::Real, Delta::Number)
    # Create the operators
    ops = Dict{String, Matrix{ComplexF32}}()

    # Empty Operators
    ops["ham"] = zeros(ComplexF32, sub.prop.dimension, sub.prop.dimension)
    ops["tb"] = zeros(ComplexF32, sub.prop.dimension, sub.prop.dimension)

    # Fock 
    st_list = fock(0, modes)
    for st in sub.prop.states
        # Position of the state
        ket = sub.prop.lookup[st]

        # Binary decomposition
        fock!(st, modes, st_list)

        # Hatugai-Kohmoto Interaction U n(k, α, ↑)n(k, α, ↓) (we use the fact that up and down spins are next to each other on the list)
        ops["ham"][ket, ket] += U * sum(st_list[2*index - 1] * st_list[2*index] for index in 1:fld(modes, 2))

        # Tight-binding operator t(α, β)c†(k, α, σ)c(k, β, σ)
        for alpha in 1:nMMHK

            # Nearest-neighbor to the right
            beta = alpha % nMMHK + 1

            for ksign in (0, 1)
                for sigma in (0, 1)
                    # Get the positions
                    nkas_pos = npos(ksign, alpha, sigma, nMMHK)
                    nkbs_pos = npos(ksign, beta, sigma, nMMHK)

                    # Get the filling
                    nkas = st_list[nkas_pos]
                    nkbs = st_list[nkbs_pos]

                    # Check if this state allows a particle to be created at alpha and destroyed at beta
                    if nkas == 0 && nkbs == 1
                        bra = sub.prop.lookup[st + 2^(modes - nkas_pos) - 2^(modes - nkbs_pos)]
                        
                        if alpha == nMMHK
                            ops["tb"][bra, ket] = (-1)^sum(st_list[(nkbs_pos + 1):(nkas_pos - 1)]; init=0) * exp(im * k * (ksign == 0 ? 1 : -1))
                        else
                            ops["tb"][bra, ket] = (-1)^sum(st_list[(nkas_pos + 1):(nkbs_pos - 1)]; init=0)
                        end
                    end
                end
            end
        end
    end

    # Add TB to the Hamiltonian (the apostrophe ' means the conjugate transpose)
    ops["ham"] += ops["tb"] + ops["tb"]'

    # Return the new subspace
    return hamSubspace(sub.prop, sub.ops, ops["ham"])
end



# ----- SOLVER -----
"""
    Solve the mean-field Hamiltonian for `L` k-points for temperature `T` and chemical potential `mu`.

    Self-consistently compute Delta starting at `delta_start` with an error of `delta_eps`.
"""
function solve(nMMHK::Int, modes::Int, L::Int, mu::Real, U::Real, g::Real, T::Real, delta_start::Real = 1e-2, delta_eps::Real = 1e-4, calpha::Real = 0.8)

    # TODO:
    # Change the kx > 0 into kx >= 0 and halve its statistical weight


    # --- GET THE K-POINTS ---

    # We have L/2 k-points in the interval (0, pi)
    kk = collect(1:fld(L, 2)) .* pi ./ (fld(L, 2))


    # --- SETUP FOCK SPACE ---

    # Generate the spaces
    fspace_prop = buildPropSubs(modes)

    # Prepare their operators
    fspace_ops = [buildOpsSubs(item, nMMHK, modes) for item in fspace_prop]

    # Debuging 
    for item in fspace_ops
        println("$(item.prop.spin) $(item.prop.states)")
        println(item.ops["nk_tk_ts"])
        println("")
    end

    # --- SOLVE ---
    # We work as follows:
    #   1. (Parallel) Solve the Hamiltonian at each k-point
    #   2. (Sync) Compute the minimum energy accross k-points, for use in thermal averages
    #   3. (Parallel) Compute expectation values for the various quantities
    #   4. (Sync) Sum the results, update self-consistent parameters and repeat 1.
    
    # We will parallelize later, the content of the parenthesis are notes to self

    # No SC term in the Hamiltonian
    if g == 0
        delta_start = 0
    end

    # Initialize the errors
    delta_error = delta_eps + 1

    # We use a flag to do a last lap after converging
    # This runs the calculation with the converged parameters
    calc_final = false
    while !calc_final

        # Check if we are converged, if we are then this is the final calculation
        if delta_error < delta_eps
            calc_final = true
        end

        # TODO: Solve the Hamiltonian
        for k in kk
            for subOps in fspace_ops
                subHam = buildHamSubs(subOps, nMMHK, modes, k, U, delta_start)
            end
        end

        delta_error = 0
    end
end

# ----- MAIN CODE -----

# Constants
W = 4

# Parameters
nMMHK = 1       # MMHK number of mixed momenta
L = 100
mu = 1.0
U = 1
g = 0
T = 0

# Auxiliary Variables
modes = 4 * nMMHK


solve(nMMHK, modes, L, mu, U, g, T)






