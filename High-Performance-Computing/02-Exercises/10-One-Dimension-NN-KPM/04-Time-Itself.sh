for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time-$trial.txt
    
    for tnum in 1 2 4 8 16 32 64 112
    do
        # Append the output to the time file
        ./a.out 100000 1 0 2 1000 2000 $tnum 0 $trial >> time-$trial.txt
    done
done
