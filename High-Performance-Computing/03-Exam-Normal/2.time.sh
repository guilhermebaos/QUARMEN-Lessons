for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time2-$trial.txt
    
    for tnum in 1 2 4 8 16
    do
        # Append the output to the time file
        ./a.out 64 90 0 0.05 0.5 600 $tnum 0 $trial >> time2-$trial.txt
    done
done
