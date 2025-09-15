for trial in 1 2 3 4
do

    echo "" > time-$trial.txt

    for tnum in 1 2 3 4 5 6 7 8
    do

        echo "Running time for tnum = $tnum" >> time-$trial.txt
        { time ./a.out 4 40000 0.005 256 128 8 1.58 0; } 2>&1 | tee -a time-$trial.txt > /dev/null
        echo " " >> time-$trial.txt

    done
done