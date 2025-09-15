for trial in 1 2 3 4
do
    # Clear the file at the start of each trial
    echo "" > time3-$trial.txt
    
    for tnum in 1 2 4 8 16
    do
        # Append the output to the time file
        ./a.out 262144 1 0.25 2 1024 $tnum 1 $trial >> time3-$trial.txt
    done
done
