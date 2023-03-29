cd src
mpic++ -c zone.cc -o temp/zone.o
mpic++ -c main.mpi.cc -o temp/main.mpi.o
mpic++ temp/zone.o temp/main.mpi.o -o nbody

#cd src
#mpic++ -c haha.cc
#mpic++ haha.o -o haha
