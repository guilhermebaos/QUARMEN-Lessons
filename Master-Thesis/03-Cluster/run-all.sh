#!/bin/bash

# Loop through every item in the current directory
for dir in */; do

    # Check if it is actually a directory
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        
        # Move into the directory
        cd "$dir" || continue

        # 1. Erase the "output" folder if it exists
        if [ -d "output" ]; then
            rm -rf "output"
        fi

        # 2. Convert line endings (useful if files were edited on Windows)
        if [ -f "run.slurm" ]; then
            dos2unix run.slurm

            # 3. Submit the job to the SLURM queue
            sbatch run.slurm
        else
            echo "  - Warning: run.slurm not found in $dir"
        fi

        # Move back up to the parent directory
        cd ..
    fi
done

echo "Done!"