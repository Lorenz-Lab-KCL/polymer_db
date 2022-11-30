#!/bin/bash -l                                                                                                                                                                                  
#$ -l h_rt=24:00:00                                                                                                                                                                             
#$ -l mem=19G                                                                                                                                                                                   
#$ -pe mpi 40                                                                                                                                                                                   
#$ -N rob_ccsd                                                                                                                                                                                  
#$ -ac allow=Y                                                                                                                                                                                  
#$ -A KCL_Admin_rse                                                                                                                                                                             
#$ -P Gold                                                                                                                                                                                      
#$ -cwd                                                                                                                                                                                         

module unload default-modules/2018
module unload compilers
module load gerun
module load gcc-libs
module load compilers/gnu/4.9.2
module load mpi/openmpi/3.1.4/gnu-4.9.2
module load orca/4.2.1-bindist/gnu-4.9.2

export OMP_NUM_THREADS=1
ORCA_EXEC=$(which orca)
$ORCA_EXEC input.inp > output.out

