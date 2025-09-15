for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time-$trial.txt
    
    for tnum in 1 2 4 8 16 32 64 100
    do
        # Append the output to the time file
        ./a.out $tnum 2000 1 256 128 4 2 0 1 $trial >> time-$trial.txt
    done
done
