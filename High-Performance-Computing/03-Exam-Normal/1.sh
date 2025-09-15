for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time1-$trial.txt
    
    for tnum in 1 2 4 8 16
    do
        # Append the output to the time file
        ./a.out $tnum 1 200000 0.05 0 $trial >> time1-$trial.txt
    done
done
