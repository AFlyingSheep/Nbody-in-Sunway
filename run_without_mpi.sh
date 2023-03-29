#bsub -I -share_size 500 -q q_share  -n 5 build/src/nbody
#bsub -I -shared -q q_share -n 2 -cgsp 64 src/haha
#bsub -I -J test -q q_share -n 2  src/haha
#bsub -I -J zhaoyang_sph_test -q q_share -n 4 -cgsp 64 -share_size 14000 src/nbody
bsub -I -J zhaoyang_sph_test -q q_share -n 1 -cgsp 64 -share_size 14000 src/nbody_without_mpi
#bsub -I -J test -q q_share -n 2 -cgsp 64 -share_size 14000 src/haha
#
