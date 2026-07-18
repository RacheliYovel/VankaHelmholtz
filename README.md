# VankaHelmholtz
Code for implementing and analyzing Vanka smoother for the acoustic Helmholtz equation, following the paper "Vanka-smoothed shifted Laplacian multigrid preconditioners for the Helmholtz equation" by Rachel Yovel, Yunhui He and Eran Treister.

The Local Fourier Analysis (LFA) code is written in Matlab and the numerical implementation in Julia language.

To run the LFA and get the optimal damping parameters for each smoother, you need to run the file "driver.m" in a folder containing the rest of the files.

To run the numerical experiments, you first need to add the following Julia packages by their git link:
https://github.com/JuliaInv/KrylovMethods.jl.git                      
https://github.com/JuliaInv/Helmholtz.jl.git               
https://github.com/JuliaInv/jInv.jl.git
(These links also appear in the documentation of the relevant file).
Then run the file "driver.jl" in a folder containing the other files. The numerical results will appear by order. For the last tables and figures, you might need to verify that your machine has at least 32 GB RAM for the code to run.
