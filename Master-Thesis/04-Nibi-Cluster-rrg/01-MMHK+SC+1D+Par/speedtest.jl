# Imports
import HDF5 as hdf

# Output a notice
println("Starting speed test!")

# Input file name
in_file = "code.jl"

# Output file name
out_file = "Cluster Test 1.test"

# Thread counts to test
tt = [1, 2, 4, 8, 16, 32, 64, 128, 192]
ss = []

# Save the data next to this file
save_path = joinpath(@__DIR__, "outputs", out_file)

# Base time, for one thread
time1 = 0.0

# Test every thread count
for t in tt
    # Launch julia with -t flag
    runtime = @elapsed run(`julia -t $t $(joinpath(@__DIR__, in_file)) 2 3600 -1 0 0 0 -1 mu -2 2 1 0.001`)
    
    if t == 1
        global time1 = runtime
    end
    
    speedup = time1 / runtime
    push!(ss, speedup)
end


# Convert arrays
ttTyped = Int64.(tt)
ssTyped = Float64.(ss)


# Save the data
hdf.h5open(save_path, "w") do file

    # Save outputs
    hdf.write(file, "tt_data", ttTyped)
    hdf.write(file, "ss_data", ssTyped)
end