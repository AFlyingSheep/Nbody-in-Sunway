cd src/master_slave
mpic++ -mftz -mieee -c main.mpe.cc -I /usr/sw/mpi/mpi_current/include/ -o build/main.mpe.o
mpic++ -mftz -mieee -c zone.cc -I /usr/sw/mpi/mpi_current/include/ -o build/zone.o
mpic++ -mftz -mieee -c calculate.cc -I /usr/sw/mpi/mpi_current/include/ -o build/calculate.o
mpic++ -mftz -mieee -c master.cc -I /usr/sw/mpi/mpi_current/include/ -o build/master.o
sw9gcc -mslave -mieee -msimd -c slave.cc -o build/slave.o

mpic++ -mhybrid -static build/calculate.o build/main.mpe.o build/zone.o build/slave.o build/master.o -o swnbody