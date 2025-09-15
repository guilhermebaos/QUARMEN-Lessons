for L in 4 8 16
do  
    for K in 25 100 225
    do
        echo "Running for L=$L and K=$K"
        ./a.out $L $K 0.0 0.02 1000 1 1 0
    done
done
