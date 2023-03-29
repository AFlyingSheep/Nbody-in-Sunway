bsub -I -J zhaoyang_sph_test -q q_share -n 20 -cgsp 64 -share_size 14000 src/nbody
#bsub -I -J zhaoyang_sph_test -q q_share -n 4 -cgsp 64 -share_size 14000 src/nbody_without_mpi
