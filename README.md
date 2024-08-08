The repository includes the codes and data files to reproduce the figures and tables in the manuscript.

Some important notes are as follows:

(1) Run_simulation.m and Run_real_data.m are the documents for reproducing the figures and tables in the manuscript.
(2) The glmnet_matlab package, which is needed when running the Run_real_data.m file, requires a Fortran compiler environment (“MinGW64 Compiler (FORTRAN)”) to be set up.
(3) Re-running the simulation program will take a considerable amount of time. If you want to directly load the results obtained from previous runs, please set “load_if = 1” in the line21 of Run_simulation.m.
(4) Please modify the address of the Transfer_fusion_regression folder on the local machine in the NewAddress at lines 20-21 of Run_simulation.m.
