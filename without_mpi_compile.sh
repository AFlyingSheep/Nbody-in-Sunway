cd src
mpic++ -c zone.cc -o temp/zone.o
mpic++ -c main.cc -o temp/main.o
mpic++ temp/zone.o temp/main.o -o nbody_without_mpi

#cd src
#mpic++ -c haha.cc
#mpic++ haha.o -o haha
