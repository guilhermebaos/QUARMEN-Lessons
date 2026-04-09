## Imports
import Base
import LinearAlgebra as la
import HDF5 as hdf

# Parallelization
using Base.Threads

### ----- HELPER FUNCTIONS -----
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
function bissect(func::Function, a::Real, b::Real, eps::Real, maxI::Int = 200, mult::Bool = false)::Real

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
    output .= ((decimal >> pos) % 2 for pos in (modes-1):-1:0)
end


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
    npos(k::Int, alpha::Int, sigma::Int, nMMHK::Int)

Returns the position of the occupation number for state (k, α, σ).

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





## ----- SUBSPACES -----
struct Subspace
    # Subspace quantum numbers
    spin::Int64

    # States in the subspace
    states::Vector{Int64}

    # Dimension of the space
    dimension::Int64

    # Lookup table for states
    lookup::Dict{Int64, Int64}

    # Operators
    # We make them all be complex because we eventually want complex multiplications
    ops::Dict{String, Matrix{ComplexF64}}

    # Eigenvalues
    vals::Vector{Float64}

    # Eigenvectors
    vecs::Matrix{ComplexF64}
end


"""
    Base.show(io::IO, sub::Subspace)

Print the properties of the subspace.
"""
function Base.show(io::IO, sub::Subspace)
    print(io, "Space with S = $(sub.spin) contains $(length(sub.states)) states: $(sub.states).")
end


"""
    startSubs(modes::Int)

Compute the elementary properties of the subspaces.
"""
function startSubs(modes::Int)
    # Dimension of the Fock space
    dim = 1 << modes

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
    fspace = [(spin, states, lookup) for (spin, states, lookup) in zip(fspace_ss, fspace_states, fspace_lookup)]

    return fspace
end






## ----- OPERATORS -----
"""
    buildOpsSubs(spin::Int, states::Vector{Int}, lookup::Dict{Int, Int}, nMMHK::Int, modes::Int)
    
Build the operators of this subspace which are independent of k and of any order parameters.
"""
function buildOpsSubs(spin::Int, states::Vector{Int}, lookup::Dict{Int, Int}, nMMHK::Int, modes::Int, mu::Real, U::Real)
    # Compute dimension
    dimension = length(states)
    
    # Create the operators
    ops = Dict{String, Matrix{ComplexF64}}()

    # Empty operators
    ops["nk_tk_ts"] = zeros(ComplexF64, dimension, dimension)

    ops["nk_tk_up"] = zeros(ComplexF64, dimension, dimension)
    ops["nk_tk_dw"] = zeros(ComplexF64, dimension, dimension)
    
    # ops["nk_pk_ts"] = zeros(Int64, dimension, dimension)
    # ops["nk_mk_ts"] = zeros(Int64, dimension, dimension)
    
    ops["bk_tk"] = zeros(ComplexF64, dimension, dimension)
    ops["fmk"] = zeros(ComplexF64, dimension, dimension)
    
    ops["hamTB"] = zeros(ComplexF64, dimension, dimension)
    ops["hamHK"] = zeros(ComplexF64, dimension, dimension)
    ops["hamMU"] = zeros(ComplexF64, dimension, dimension)


    # Fill in the operators
    st_list = fock(0, modes)
    for st in states
        # Position of the state
        ket = lookup[st]

        # Binary decomposition, by mutating st_list
        fock!(st, modes, st_list)

        # Number of particles
        ops["nk_tk_ts"][ket, ket] = sum(st_list)

        # Number of particles with a given spin (odds are spin up, evens are spin down)
        ops["nk_tk_up"][ket, ket] = sum(n * (pos % 2) for (pos, n) in enumerate(st_list))
        ops["nk_tk_dw"][ket, ket] = sum(n * ((pos + 1) % 2) for (pos, n) in enumerate(st_list))

        # Number of particles at +k and -k
        # ops["nk_pk_ts"][ket, ket] = sum(n * (pos <= modes / 2) for (pos, n) in enumerate(st_list))
        # ops["nk_mk_ts"][ket, ket] = sum(n * (pos >  modes / 2) for (pos, n) in enumerate(st_list))


        # Chemical Potential H = -μN
        ops["hamMU"][ket, ket] = -mu * sum(st_list)

        # Hatugai-Kohmoto Interaction H = U n(k, α, ↑)n(k, α, ↓) (we use the fact that up and down spins are next to each other on the list)
        ops["hamHK"][ket, ket] = U * sum(st_list[2*index - 1] * st_list[2*index] for index in 1:(2*nMMHK))


        # Applying bk to this state for each alpha (st is the ket, we find the bra)
        # The operator is given by bk = c(-p, α, ↓)c(p, α, ↑)
        for alpha in 1:nMMHK
            # Applying bk for +k
            # if nk_mk_ds * nk_pk_ds for this orbital
            pos_pk_us = npos(0, alpha, 0, nMMHK)
            pos_mk_ds = npos(1, alpha, 1, nMMHK)
            if (st_list[pos_pk_us] * st_list[pos_mk_ds]) != 0
                # The bra is the state without those two particles
                bra_pk = lookup[st - 1 << (modes - pos_pk_us) - 1 << (modes - pos_mk_ds)]

                # The phase is given by the states between the destroyed particles
                # We use iseven(x) instead of (-1)^x
                ops["bk_tk"][bra_pk, ket] = iseven(sum(st_list[(pos_pk_us + 1):(pos_mk_ds - 1)]; init=0)) ? 1 : -1
            end

            # Applying bk for -k
            # if nk_pk_ds * nk_mk_ds for this orbital
            pos_pk_ds = npos(0, alpha, 1, nMMHK)
            pos_mk_us = npos(1, alpha, 0, nMMHK)
            if (st_list[pos_pk_ds] * st_list[pos_mk_us]) != 0
                # The bra is the state without those two particles (the -1 is because Julia is 1-indexed)
                bra_mk = lookup[st - 1 << (modes - pos_pk_ds) - 1 << (modes - pos_mk_us)]

                # The phase is given by the states between the destroyed particles
                # The +1 is from the fact that now we destroy the -k first, which has to commute with the filled +k state
                ops["bk_tk"][bra_mk, ket] = iseven(sum(st_list[(pos_pk_ds + 1):(pos_mk_us - 1)]; init=0) + 1) ? 1 : -1
            end
        end
    end

    # Magnetization
    # TODO: Make magnetization per orbital
    ops["fmk"] .= ops["nk_tk_up"] - ops["nk_tk_dw"]

    # Create empty arrays for vals and vecs
    vals = Vector{Float64}(undef, dimension)
    vecs = Matrix{ComplexF64}(undef, dimension, dimension)

    # Return the subspace
    return Subspace(spin, states, dimension, lookup, ops, vals, vecs)
end


"""
    buildHamTB!(sub::Subspace, nMMHK::Int, modes::Int, k::Real)

Create the tight-binding term of the Hamiltonian.
"""
function buildHamTB!(sub::Subspace, nMMHK::Int, modes::Int, k::Real)

    # In this case the tight-binding term is just the usual dispersion
    if nMMHK == 1
        sub.ops["hamTB"] .= ComplexF64.(2 * cos(k) .* sub.ops["nk_tk_ts"])
        return
    end

    # Empty Operator
    sub.ops["hamTB"] .= ComplexF64.(0)

    # Fock 
    st_list = fock(0, modes)
    for st in sub.states
        # Position of the state
        ket = sub.lookup[st]

        # Binary decomposition
        fock!(st, modes, st_list)

        # Tight-binding operator H = t(α, β)c†(k, α, σ)c(k, β, σ)
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
                        bra = sub.lookup[st + 1 << (modes - nkas_pos) - 1 << (modes - nkbs_pos)]
                        
                        if alpha == nMMHK
                            sub.ops["hamTB"][bra, ket] += (iseven(sum(st_list[(nkbs_pos + 1):(nkas_pos - 1)]; init=0)) ? 1 : -1)  * exp(im * k * (ksign == 0 ? 1 : -1))
                        else
                            sub.ops["hamTB"][bra, ket] += iseven(sum(st_list[(nkas_pos + 1):(nkbs_pos - 1)]; init=0)) ? 1 : -1
                        end
                    end
                end
            end
        end
    end

    # Add the adjoint TB to the Hamiltonian (the apostrophe ' means the conjugate transpose)
    sub.ops["hamTB"] .= sub.ops["hamTB"] .+ sub.ops["hamTB"]'
end





## ----- OBSERVABLES -----
"""
    thermal_average(fspace::Array{Subspace}, opcode::String, T::Real, make_positive::Bool = false)

Compute the thermal average of an operator `opcode`.

# Arguments
- `fspace::Array{Subspace}`: The list with the subspaces.
- `vecs::Array{Number}`: Array with `2^(modes)` eigenvectors.
- `vals::Array{Number}`: Array with `2^(modes)` eigenvalues.

The arrays `vecs` and `vals` are ordered in blocks of `fspace[i].dimension` elements, each relative to their `fspace[i]`.
"""
function thermal_average(fspace::Array{Subspace}, opcode::String, T::Real, vals::Vector{Float64}, make_positive::Bool = false)

    # Join all the eigenvalues
    index = 1
    for sub in fspace
        newIndex = index + sub.dimension - 1
        vals[index:newIndex] = sub.vals
        index = newIndex + 1
    end

    # Scaled energies
    mini = minimum(vals)

    # Boltzmann exponential
    # If T = 0 then if E = 0 the exponent is 1 and otherwise is zero
    if isapprox(T, 0.0, atol=1e-14)
        for i in eachindex(vals)
            vals[i] = isapprox((vals[i] - mini), 0.0, atol=1e-14) ? 1.0 : 0.0
        end
    else
        for i in eachindex(vals)
            vals[i] = exp(-(vals[i] - mini) / T)
        end
    end
    bbexp = vals
    
    # Partition function
    Z = sum(bbexp)

    # Boltzmann exponential -> Boltzmann factor
    bbexp ./= Z

    # Compute average via a trace
    result = 0
    vecnum = 1
    for sub in fspace
        
        # Choose operator
        op = sub.ops[opcode]

        # Average = < n | O | n > = O(i, j) n*(i) n(j)
        # Then we multiply it by the weight exp(-En / T)
        # Finally we take the trace by summing over all eigenstates n
        for n in 1:sub.dimension 
            average = 0
            for j in 1:sub.dimension
                v_jn = sub.vecs[j, n]

                for i in 1:sub.dimension
                    # Debugging
                    # println("$n, $i, $j $(vecs[n])")

                    # Compute the average
                    average += op[i, j] * sub.vecs[i, n]' * v_jn
                end
            end
            
            # Compute the absolute value of the expectation value on each state
            if make_positive
                average = abs(average)
            end

            result += average * bbexp[vecnum]
            vecnum += 1
        end
    end

    return result
end





## ----- SOLVER -----
"""
    solve(nMMHK::Int, L::Int, mu::Real, U::Real, g::Real, T::Real, delta_start::Real = 1e-2, delta_eps::Real = 1e-4, calpha::Real = 0.8)

Solve the mean-field Hamiltonian for `L` k-points for temperature `T` and chemical potential `mu`.

Self-consistently compute Delta starting at `delta_start` with an error of `delta_eps`.
"""
function solve(nMMHK::Int, L::Int, mu::Real, U::Real, g::Real, T::Real, delta_start::ComplexF64 = 1e-2 + 0 * im, delta_eps::Real = 1e-4, calpha::Real = 0.8)

    # TODO:
    # Change the kx > 0 into kx >= 0 and halve its statistical weight
    
    # --- Auxiliary variables ---
    modes = 4 * nMMHK



    # --- GET THE K-POINTS ---

    # We have L/(2 * nMMHK) k-points in the interval (0, pi / nMMHK)
    kk = collect(1:fld(L, 2 * nMMHK)) .* pi ./ (fld(L, 2 * nMMHK))

    # Compute actual number of k-points
    Nk = length(kk)


    # --- SETUP FOCK SPACE ---

    # Generate the spaces
    fspace = startSubs(modes)

    # Prepare their operators
    fspace = [buildOpsSubs(spin, states, lookup, nMMHK, modes, mu, U) for (spin, states, lookup) in fspace]

    # Debuging 
    # for item in fspace
    #     println("$(item.spin) $(item.states)")
    #     println(item.ops["nk_tk_ts"])
    #     println("")
    # end



    # Initialize the errors
    delta_error = delta_eps + 1

    # No SC term in the Hamiltonian
    if g == 0
        delta_start = 0
        delta_error = delta_eps - 1
    end



    # --- SOLVE ---
    # We work as follows:
    #   1. (Parallel) Solve the Hamiltonian at each k-point
    #   2. (Parallel) Compute the minimum energy at each k-points, for use in thermal averages
    #   3. (Parallel) Compute expectation values for the various quantities
    #   4. (Sync) Sum the results, update self-consistent parameters and repeat 1.

    # Outputs
    n_atomic = Atomic{Float64}(0)
    Kxx_atomic = Atomic{Float64}(0)

    # Split the k-points in chunks, one for each thread
    chunk_size = cld(Nk, Threads.nthreads())
    chunks = collect(Iterators.partition(kk, chunk_size))

    # Compute until the error is smaller then the desired precision
    keep_going = true
    while keep_going
        
        # This is the last lap!
        if delta_error < delta_eps
            keep_going = false
        end

        # Multi-thread variable for the order paramter
        # Atomic opeartions are not defined for complex numbers, so we separate into real and imaginary parts
        delta_atomic_real = Atomic{Float64}(0)
        delta_atomic_imag = Atomic{Float64}(0)

        # [PAR] Start paralelization (each thread needs a deepcopy of fspace)
        @threads for ch in chunks

            # Copiying the needed data
            fspace_local = deepcopy(fspace)
            valsOverwrite = Vector{Float64}(undef, 2^modes)
            
            # Accumulators
            deltaLocal = 0.0 + 0.0im
            nLocal = 0.0
            KxxLocal = 0.0
        
            for k in ch

                # Solve the system in each subspace
                for sub in fspace_local

                    # Build the Hamiltonian for this k in each subspace
                    buildHamTB!(sub, nMMHK, modes, k)

                    # Overwrite the tight-binding matrix for efficiency
                    sub.ops["hamTB"] .+= sub.ops["hamMU"] .+ sub.ops["hamHK"] .+ (delta_start' .* sub.ops["bk_tk"] .+ delta_start .* sub.ops["bk_tk"]')

                    # Solve it
                    eigensolved = la.eigen!(la.Hermitian(sub.ops["hamTB"]))
                    sub.vals .= eigensolved.values
                    sub.vecs .= eigensolved.vectors
                end

                # Compute Delta
                deltaLocal += thermal_average(fspace_local, "bk_tk", T, valsOverwrite)

                # This is the last lap, compute outputs
                if !keep_going
                    nLocal += real(thermal_average(fspace_local, "nk_tk_ts", T, valsOverwrite))

                    # Get the pure kinetic Hamiltonian
                    for sub in fspace_local
                        buildHamTB!(sub, nMMHK, modes, k)
                    end

                    KxxLocal += real(thermal_average(fspace_local, "hamTB", T, valsOverwrite))
                end
            end

            # [ONE-THREAD] Add the local value to the global total
            atomic_add!(delta_atomic_real, real(deltaLocal))
            atomic_add!(delta_atomic_imag, imag(deltaLocal))
            if !keep_going
                atomic_add!(n_atomic, nLocal)
                atomic_add!(Kxx_atomic, KxxLocal)
            end
        end

        # Get complex delta_new
        delta_new = delta_atomic_real[] + delta_atomic_imag[] * im

        # Compute error and next iteration
        delta_new *= -g / (2 * Nk * nMMHK)
        delta_error = abs(delta_new - delta_start)
        delta_start = calpha * delta_new + (1 - calpha) * delta_start
    end

    # Make the outputs non-atomic
    n = n_atomic[]
    Kxx = Kxx_atomic[]

    # Normalize
    n *= 1 / (2 * Nk * nMMHK)
    Kxx *= -pi / (2 * Nk * nMMHK)

    outputs = Dict{String, Number}(
        "Nk" => Nk, 
        "n" => n,
        "Kxx" => Kxx,
        "Delta" => delta_start
    )

    return outputs
end



## ----- RUN THE SOLVER -----
function compute(params::Dict{String, Real})

    # Go for target filling
    if params["nTarget"] >= 0
        # Get minimum and maximum from the parameters
        mu_min = params["mu_min"]
        mu_max = params["mu_max"]
        mu_eps = params["mu_eps"]

        # Bissect for mu
        mu = bissect(mu_test -> solve(params["nMMHK"], params["L"], mu_test, params["U"], params["g"], params["T"])["n"] - params["nTarget"], mu_min, mu_max, mu_eps)

        # Update mu
        params["mu"] = mu
    end

    # File title
    file_arr = []
    for (key, value) in params
        if key in ["runID", "mu", "mu_min", "mu_max", "mu_eps", "nTarget"]
            continue
        else
            push!(file_arr, "$key=$value")
        end
    end
    file_str = string(params["runID"]) * "-" * join(file_arr, "-") * ".h5"

    # Compute and get the output
    output = solve(params["nMMHK"], params["L"], params["mu"], params["U"], params["g"], params["T"])

    # Save the data into this file
    save_path = joinpath(@__DIR__, "outputs", file_str)

    # Save the data
    hdf.h5open(save_path, "w") do file
        
        # Save the parameters
        g_params = hdf.create_group(file, "params")
        for (k, v) in params
            hdf.write(g_params, k, v)
        end
        
        # Save the outputs
        g_out = hdf.create_group(file, "out")
        for (k, v) in output
            hdf.write(g_out, k, v)
        end

    end
end



## ----- MAIN CODE -----
runID = parse(Int64, ARGS[1])
nMMHK = parse(Int64, ARGS[2])
L = parse(Int64, ARGS[3])
mu = parse(Float64, ARGS[4])
U = parse(Float64, ARGS[5])
g = parse(Float64, ARGS[6])
T = parse(Float64, ARGS[7])
nTarget = parse(Float64, ARGS[8])

params = Dict{String, Real}(
    "runID"   => runID,
    "nMMHK"   => nMMHK,
    "L"       => L,
    "mu"      => mu,
    "U"       => U,
    "g"       => g,
    "T"       => T,
    "nTarget" => nTarget,
    "mu_min"  => -2.2 - T,
    "mu_max"  => +2.2 + U + g + T,
    "mu_eps"  => 0.0001,
)

@time compute(params)
