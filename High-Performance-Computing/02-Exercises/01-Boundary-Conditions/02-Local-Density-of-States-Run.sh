for N in 64 256 1024
do

    ./a.out 4 0 $N 0.01 1 1
    ./a.out 4 1 $N 0.01 1 1

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