## --- IMPORTS ---
using LinearAlgebra
using Plots
using Plots.Measures


## --- CONVENTIONS ---
# We use a = 1


## --- AUXILIARY FUNCTION ---

"""
Compute the function `f` using `a = 1`.
"""
function func(kx::Real, ky::Real)
    return 1 + exp(-im * (sqrt(3) * kx + 3 * ky) / 2) + exp(-im * (-sqrt(3) * kx + 3 * ky) / 2)
end


function delxfunc(kx::Real, ky::Real)
    return -im * sqrt(3) / 2 * exp(-im * (sqrt(3) * kx + 3 * ky) / 2) + im * sqrt(3) / 2 * exp(-im * (-sqrt(3) * kx + 3 * ky) / 2)
end


function delyfunc(kx::Real, ky::Real)
    return -im * 3 / 2 * exp(-im * (sqrt(3) * kx + 3 * ky) / 2) - im * 3 / 2 * exp(-im * (-sqrt(3) * kx + 3 * ky) / 2)
end




## --- HAMILTONIAN ---

"""
Create the minimal tight-binding Hamiltonian for n-layer rhombohedral graphene.

# Arguments
- 

### Returns
- `hamArray`: Overwriten matrix hamArray with the new Hamiltonian in the upper-triangular part.
"""
function hamMin!(n::Int64, gam0::Real, gam1::Real, V::Real, kx:: Real, ky::Real, hamArray::Matrix{Complex{Float64}})
    # Write the Hamiltonian into hamArray
    hamArray .= 0
    for layer in 1:n
        # Relevant sites
        Asite = 2*layer - 1
        Bsite = 2*layer

        # Local energy, from the displacement field
        # Each layer has a difference in energy V from its neighbors
        hamArray[Asite, Asite] = hamArray[Bsite, Bsite] = (layer - (n + 1) / 2) * V

        # Intra-layer hopping γ0
        hamArray[Asite, Bsite] = -gam0 * func(kx, ky)

        # Inter-layer hopping γ1
        if layer != n
            hamArray[Bsite, Asite + 2] = gam1
        end
    end
end



"""
Create the minimal tight-binding Hamiltonian for n-layer rhombohedral graphene, differentiated along kx.


### Returns
- `hamArray`: Overwriten matrix hamArray with the new Hamiltonian (not just upper-triangular!).
"""
function delxhamMin!(n::Int64, gam0::Real, kx:: Real, ky::Real, hamArray::Matrix{Complex{Float64}})
    # Write the Hamiltonian into hamArray
    hamArray .= 0
    for layer in 1:n
        # Relevant sites
        Asite = 2*layer - 1
        Bsite = 2*layer

        # Intra-layer hopping γ0
        hamArray[Asite, Bsite] = -gam0 * delxfunc(kx, ky)
        hamArray[Bsite, Asite] = -gam0 * conj(delxfunc(kx, ky))
    end
end


function delyhamMin!(n::Int64, gam0::Real, kx:: Real, ky::Real, hamArray::Matrix{Complex{Float64}})
    # Write the Hamiltonian into hamArray
    hamArray .= 0
    for layer in 1:n
        # Relevant sites
        Asite = 2*layer - 1
        Bsite = 2*layer

        # Intra-layer hopping γ0
        hamArray[Asite, Bsite] = -gam0 * delyfunc(kx, ky)
        hamArray[Bsite, Asite] = -gam0 * conj(delyfunc(kx, ky))
    end
end



function solveMin!(n::Int64, gam0::Real, gam1::Real, V::Real, kkx::Vector{Float64}, kky::Vector{Float64}, outEnergies::Matrix{Float64}, omega::Vector{Float64}, bandN::Int64 = -1)

    # Create memory for the Hamiltonian
    ham = zeros(ComplexF64, 2*n, 2*n)

    # Create memory for the Hamiltonians
    delxham = zeros(ComplexF64, 2*n, 2*n)
    delyham = zeros(ComplexF64, 2*n, 2*n)


    for index in eachindex(kkx)

        # Get the k-point
        kx = kkx[index]
        ky = kky[index]

        # Build the Hamiltonian at this k-point
        hamMin!(n, gam0, gam1, V, kx, ky, ham)

        # Solve the Hamiltonian
        # TODO: Use mutation here to decrease allocations
        eigensolved = eigen!(Hermitian(ham))

        # Write the eigenvalues to the output
        # They are in increasing order, so we should get the correct band index
        outEnergies[:, index] .= eigensolved.values
        vecs = eigensolved.vectors

        # Compute Berry curvature
        if bandN > 0

            # Compute the Hamiltonians
            delxhamMin!(n, gam0, kx, ky, delxham)
            delyhamMin!(n, gam0, kx, ky, delyham)

            # View of the eigenvector for band n
            vn = @view vecs[:, bandN]

            omega[index] = 0.0
            for bandM in 1:2n
                if bandM != bandN
                    # View of the eigenvector for band m
                    vm = @view vecs[:, bandM]

                    # Do the < n | ∂xH | m >
                    tempx = vn' * delxham * vm

                    # Do the < m | ∂yH | n >
                    tempy = vm' * delyham * vn

                    omega[index] += -2 * imag(tempx * tempy) / (eigensolved.values[bandM] - eigensolved.values[bandN])^2
                end
            end
        end
    end
end



## --- COMPUTE ---
function plotBands(L::Int64, n::Int64, gam0::Real, gam1::Real, V::Real, viewK::Real, limy::Real, disp::Bool = false)
    
    # Dirac Point
    Diracx = 4 * sqrt(3) * pi / 9 

    # Get bands along the ky = 0 direction 
    bands = zeros(Float64, 2*n, L)
    kkx = zeros(Float64, L) 
    kky = zeros(Float64, L) 
    omega = zeros(Float64, L) 

    # Linear space in kx
    kkx = collect(0:L-1) / L * 4 * pi * sqrt(3) / 3

    # Solve for the parameters
    @time solveMin!(n, gam0, gam1, V, kkx, kky, bands, omega)

    # Initialize the plot
    plt = plot(xlabel = "k_x / K", ylabel = "Energy", title = "\nBands for ky = 0 near K\n With n = $n, V = $V\n\n")

    # Loop through each band to plot and label
    for i in 1:2n
        plot!(plt, kkx / Diracx, bands[i, :], linewidth=2, label="")
    end

    if disp
        display(plt)
    end

    # Zoom around the Dirac point
    xlims!(1 - 2*viewK, 1 + 2*viewK)
    ylims!(-limy*3, +limy*3)
    if disp
        display(plt)
    end
    
    # Zoom even more!
    Diracx = 4 * sqrt(3) * pi / 9 
    xlims!(1 - viewK, 1 + viewK)
    ylims!(-limy, +limy)
    if disp
        display(plt)
    end

    return plt
end


function plotBerry(L::Int64, n::Int64, gam0::Real, gam1::Real, V::Real, viewK::Real, bandN::Int64)

    # Dirac Point
    Diracx = 4 * sqrt(3) * pi / 9 

    # Get uniform k points in a square
    kkx = (collect(0:L) ./ L .- 0.5) .* 2viewK .* Diracx # .+ Diracx
    kky = (collect(0:L) ./ L .- 0.5) .* 2viewK .* Diracx

    gridkx = vec([x for x in kkx, y in kky])
    gridky = vec([y for x in kkx, y in kky])

    # Prepare memory
    bands = zeros(Float64, 2*n, (L+1)^2)
    omega = zeros(Float64, (L+1)^2)

    # Solve for the parameters
    @time solveMin!(n, gam0, gam1, V, gridkx, gridky, bands, omega, bandN)

    # Reshape omega into a matrix
    omega_matrix = reshape(omega, length(kky), length(kkx))

    # Divide by the maximum value
    omega_matrix /= maximum(abs.(omega_matrix))

    # Initialize the plot
    return heatmap(kkx / Diracx, kky / Diracx, omega_matrix, xlabel="kx / K", ylabel="ky / K", colorbar_title = "Ωn / Ωmax", colorbar_width = 3, title = "\nBerry of band $bandN \n With n = $n, V = $V\n\n")
end




## --- Parameters ---
L = 350
n = 5
gam0 = 2600
gam1 = 360
V = 25
viewK = 1.1
bandN = n + 1
limy = 500



# Plot many Berry curvatures and bands
hh = []
pp = []

nn = collect(3:1:7)
for paramN in nn
    hm = plotBerry(L, paramN, gam0, gam1, V, viewK, paramN + 1)
    push!(hh, hm)

    pl = plotBands(20*L, paramN, gam0, gam1, V, viewK, limy)
    push!(pp, pl)
end

plot(pp..., layout = (1, length(nn)), size=(2000, 500), left_margin = 5mm, right_margin = 5mm, top_margin = 10mm, bottom_margin = 10mm)
plot(hh..., layout = (2, div(length(nn), 2) + 1), size=(2000, 900), left_margin = 5mm, right_margin = 5mm, top_margin = 10mm, bottom_margin = 10mm)




# Plot many Berry curvatures and bands
hh = []
pp = []

vv = collect(0:25:200)
for paramV in vv
    hm = plotBerry(L, n, gam0, gam1, paramV, viewK, bandN)
    push!(hh, hm)

    pl = plotBands(20*L, n, gam0, gam1, paramV, viewK, limy)
    push!(pp, pl)
end

plot(pp..., layout = (div(length(vv), 4) + 1, 4), size=(2200, 300 * (div(length(vv), 4) + 1)), left_margin = 5mm, right_margin = 5mm, top_margin = 8mm)
plot(hh..., layout = (div(length(vv), 4) + 1, 4), size=(2200, 300 * (div(length(vv), 4) + 1)), left_margin = 5mm, right_margin = 5mm, top_margin = 8mm)




## --- OLD CODE ---

# function solveMin!(L::Int64, n::Int64, gam0::Real, gam1::Real, V::Real, outEnergies::Array{Float64, 3}, outKx::Matrix{Float64}, outKy::Matrix{Float64})
#     # Reciprocal lattice vectors
#     b1 = 4 * pi / 3 * [sqrt(3)/2, 1/2]
#     b2 = 4 * pi / 3 * [-sqrt(3)/2, 1/2] * (-1)

#     # K-vector indices
#     m1 = collect(0:L - 1) / L
#     m2 = collect(0:L - 1) / L

#     # Create memory for the Hamiltonian
#     ham = zeros(ComplexF64, 2*n, 2*n)
    
#     for (index1, coef1) in enumerate(m1)
#         for (index2, coef2) in enumerate(m2)

#             # Get the k-point
#             kx = b1[1] * coef1 + b2[1] * coef2
#             ky = b1[2] * coef1 + b2[2] * coef2

#             # Save the k-point
#             outKx[index1, index2] = kx
#             outKy[index1, index2] = ky

#             # Build the Hamiltonian at this k-point
#             hamMin!(n, gam0, gam1, V, kx, ky, ham)

#             # Solve the Hamiltonian
#             eigensolved = eigen!(Hermitian(ham))

#             # Write the eigenvalues to the output
#             # They are in increasing order, so we should get the correct band index
#             outEnergies[:, index1, index2] .= eigensolved.values
#         end
#     end
# end


# function plotBands(L::Int64, n::Int64, gam0::Real, gam1::Real, V::Real, dirac::Real, limy)
    
#     # Create empty arrays
#     energies = zeros(Float64, 2*n, L, L)
#     kkx = zeros(Float64, L, L)
#     kky = zeros(Float64, L, L)

#     # Solve for the parameters
#     @time solveMin!(L, n, gam0, gam1, V, energies, kkx, kky)

#     # Get bands along the ky = 0 direction 
#     bands = zeros(Float64, 2*n, L)
#     kk = zeros(Float64, L) 

#     k_index = 1
#     for m1 in 1:L
#         for m2 in 1:L
#             if isapprox(kky[m1, m2], 0, atol=1e-10)
#                 bands[:, k_index] = energies[:, m1, m2]
#                 kk[k_index] = kkx[m1, m2]
#                 k_index += 1
#             end
#         end
#     end
#     # Sort the data to get good lines
#     mask = sortperm(kk)
#     kk_sorted = kk[mask]
#     bands_sorted = bands[:, mask]


#     # Initialize the plot
#     plt = plot(xlabel = "k_x", ylabel = "Energy", title = "Bands (ky = 0, γ0 = $gam0, γ1 = $gam1, V = $V)", legend = :outerright)

#     # Loop through each band to plot and label
#     for i in 1:2n
#         plot!(plt, kk_sorted, bands_sorted[i, :], linewidth=2)
#     end
#     display(plt)

#     # Zoom around the Dirac point
#     Diracx = 4 * sqrt(3) * pi / 9 
#     xlims!(Diracx - 2*dirac, Diracx + 2*dirac)
#     ylims!(-limy*3, +limy*3)
#     display(plt)
    
#     # Zoom even more!
#     Diracx = 4 * sqrt(3) * pi / 9 
#     xlims!(Diracx - dirac, Diracx + dirac)
#     ylims!(-limy, +limy)
#     display(plt)
# end