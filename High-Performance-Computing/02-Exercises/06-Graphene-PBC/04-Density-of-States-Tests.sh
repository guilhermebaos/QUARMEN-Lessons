for N in 17 34
do
    for gam in 0.05 0.15 
    do
        ./a.out 1 $N $gam 0 0 1 3.5 1000 1;
    done
done