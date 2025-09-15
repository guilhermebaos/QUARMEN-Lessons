for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time-$trial.txt
    
    for tnum in 1 4 9 16 25 36 49 64 81 100
    do
        # Append the output to the time file
        ./a.out 300 1 0 2 1024 $tnum 0 $trial >> time-$trial.txt
    done
done
