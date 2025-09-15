echo "" > time.txt

for tnum in 1 2 3 4 5 6 7 8
do

    echo "Running time for tnum = $tnum" >> time.txt
    { time ./a.out $tnum 1 256 0.01 0; } 2>&1 | tee -a time.txt
    echo " " >> time.txt

done


# for tnum in 1 2 4
# do

#     for pbc in 0 1
#     do

#         for N in 64 128 4096 8192
#         do
            
#             for gam in 0.1 0.01 0.001
#             do

#                 ./a.out 4 1 $N 0.01

#             done
#         done
#     done
# done