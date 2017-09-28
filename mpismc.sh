#!/bin/bash

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml NULL res2_mpi_smc_14.txt 1000 14

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_14.txt res2_mpi_smc_21.txt 1000 21

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_21.txt res2_mpi_smc_28.txt 1000 28

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_28.txt res2_mpi_smc_35.txt 1000 35

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_35.txt res2_mpi_smc_42.txt 1000 42

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_42.txt res2_mpi_smc_49.txt 1000 49

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_49.txt res2_mpi_smc_56.txt 1000 56

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_56.txt res2_mpi_smc_63.txt 1000 63

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_63.txt res2_mpi_smc_70.txt 1000 70

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_70.txt res2_mpi_smc_77.txt 1000 77

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_77.txt res2_mpi_smc_84.txt 1000 84

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_84.txt res2_mpi_smc_91.txt 1000 91

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_91.txt res2_mpi_smc_105.txt 1000 105

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_105.txt res2_mpi_smc_119.txt 1000 119

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_119.txt res2_mpi_smc_133.txt 1000 133

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_133.txt res2_mpi_smc_147.txt 1000 147

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_147.txt res2_mpi_smc_161.txt 1000 161

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_161.txt res2_mpi_smc_175.txt 1000 175

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_175.txt res2_mpi_smc_189.txt 1000 189

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_189.txt res2_mpi_smc_203.txt 1000 203

mpirun -np 8 ./mpismc/mpismc ./mpiabck/problem.xml res2_mpi_smc_203.txt res2_mpi_smc_217.txt 1000 217