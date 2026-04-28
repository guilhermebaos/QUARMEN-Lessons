using HDF5

function aggregate_hdf5_files(input_dir::String, output_file::String)
    # Get all .h5 files in the target directory
    files = filter(f -> endswith(f, ".h5"), readdir(input_dir, join=true))
    
    # Check that some files were found
    if isempty(files)
        println("No HDF5 files found in the directory: $input_dir")
        return
    end
    println("Found $(length(files)) files. Starting aggregation...")

    # Open the new aggregated file in write mode
    h5open(output_file, "w") do h5_out
        
        for file in files

            # The filename without the path and extension is the unique group name
            run_group_name = splitext(basename(file))[1]
            println("  -> Aggregating: $run_group_name")
            
            # Open the individual run file in read mode
            h5open(file, "r") do h5_in

                # Create a top-level group for this specific run in the master file
                g_run = create_group(h5_out, run_group_name)
                
                # Iterate over the top-level groups in the input file (should be "params", "out")
                for group_name in keys(h5_in)
                    g_src = h5_in[group_name]
                    g_dst = create_group(g_run, group_name)
                    
                    # Read each dataset/attribute and write it to the new file
                    for dataset_name in keys(g_src)
                        data = read(g_src[dataset_name])
                        write(g_dst, dataset_name, data)
                    end
                end
            end
        end
    end
    
    println("\nAggregation complete! All data saved to: $output_file")
end

# Define paths relative to the script's location
input_directory = joinpath(@__DIR__, "outputs")
output_filepath = joinpath(@__DIR__, "output-agg.h5")

# Run the aggregation
aggregate_hdf5_files(input_directory, output_filepath)