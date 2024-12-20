# Axonal-PathFinding
This repository contains the ABAQUS VUMAT subroutine code for modeling axonal growth and pathfinding during brain development. The folder "Abaqus_inp_sample" includes a sample geometrical input file. This code has been tested on Abaqus 2019. To run the Abaqus job, use the following command:

Abaqus job=G_6_Case_0.inp user=vumat_axongrowth_G_6.f cpus=32 thread_mode=MPI interactive

You may adjust the number of CPUs as needed. Expected run time for a model with 100 fibers on a "normal" desktop computer with 32 cpus should be around two hours. Before launching the job, ensure that Fortran is properly linked with ABAQUS. Refer to https://gist.github.com/franaudo/72362784ded685e4cb381e57020c9ec7 for instructions. Additionally, this code utilizes MPI to accelerate the running time, so the MPI library should also be installed on your system. For installation instructions, please refer to https://github.com/mwierszycki/openmpi-abaqus/)

If you use our code in your research or project, please cite it by referencing the DOI 10.5281/zenodo.14536988 assigned to this repository. Proper citation helps us track the impact of our work and acknowledges the effort that went into developing this code.
