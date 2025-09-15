for L in 4 16 64
do
    for w in 0.25 0.5 1 2 4
    do
        ./a.out $L 120 $w 0.05 0.5 600 16 1 0
    done
done