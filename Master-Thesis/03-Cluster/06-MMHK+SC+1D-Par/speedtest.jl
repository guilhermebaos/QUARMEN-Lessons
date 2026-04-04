# Imports
import HDF5 as hdf

# Input file name
in_file = "code.jl"

# Output file name
out_file = "nMMHK = 2 L = 30000 muSte = 4.test"

# Thread counts to test
tt = [1, 2, 4]
ss = []

# Create outputs folder
mkpath(joinpath(@__DIR__, "outputs"))

# Save the data next to this file
save_path = joinpath(@__DIR__, "outputs", out_file)

# Base time, for one thread
time1 = 0.0

# Test every thread count
for t in tt
    # Launch julia with -t flag
    runtime = @elapsed run(`julia -t $t $(joinpath(@__DIR__, in_file))`)
    
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