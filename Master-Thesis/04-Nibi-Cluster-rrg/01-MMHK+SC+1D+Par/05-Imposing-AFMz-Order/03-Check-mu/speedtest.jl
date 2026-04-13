# Imports
import HDF5 as hdf

# Input file name
in_file = "code.jl"

# Output file name
out_file = "speed.test"

# Thread counts to test
tt = [32, 16, 8, 4, 2, 1]
rr = []

# Save the data next to this file
save_path = joinpath(@__DIR__, out_file)

# Base time, for one thread
time1 = 0.0

# Test every thread count
for t in tt
    # Launch julia with -t flag
    runtime = @elapsed run(`julia -t $t $(joinpath(@__DIR__, in_file)) $t 2 1920000 0 0 0 0 -1`)
    push!(rr, runtime)
end


# Convert arrays
ttTyped = Int64.(tt)
rrTyped = Float64.(rr)


# Save the data
hdf.h5open(save_path, "w") do file

    # Save outputs
    hdf.write(file, "tt_data", ttTyped)
    hdf.write(file, "rr_data", rrTyped)
end