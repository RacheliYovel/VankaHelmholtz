# VankaHelmholtz
Code for reproducing the experiments in the paper "Vanka-smoothed shifted Laplacian multigrid preconditioners for the Helmholtz equation" by Rachel Yovel, Yunhui He and Eran Treister.

The code contains two main folders: 
- Thoretical local Fourier analysis (LFA) written in Matlab,
- Numerical experiments written in Julia.

The numerical experiments depend on several formal Julia packages (can be added directly from Julia), and the following packages available in the links:
https://github.com/JuliaInv/KrylovMethods.jl.git                      
https://github.com/JuliaInv/Helmholtz.jl.git               
https://github.com/JuliaInv/jInv.jl.git    
In the files 
Manifest.toml
and
Project.toml
you can find a pinned Juia environment for which the code should run.

# Instructions for the reproduction of each result in the paper

- Table 1: run the file Numerical - Julia\experiments\table1.jl and the resulting table would appear as Numerical - Julia\output\Table1.csv

- Table 2 consists of hyperparameters, most achieved by trial and error (and some overlap with the results from Figure 5) and hence it is not reproduced by the code.

- Table 3: run the file Numerical - Julia\experiments\table3.jl and the resulting table would appear as Numerical - Julia\output\Table3.csv

- Table 4: run the file Numerical - Julia\experiments\table4.jl and the resulting table would appear as Numerical - Juliaoutput\Table4.csv

- Table 5: run the file Numerical - Julia\experiments\table5.jl and the resulting table would appear as Numerical - Julia\output\Table5.csv
* note that this table contains timing results that can very between dofferent machines and experiments

- Table 6: This experiment require using the overthrust model (SEG), first read the instructions on the end of this document and run the experiment only after the file Numerical - Julia\experiments\BenchmarkModels\3DOverthrust801801187.dat is saved.
Next, run the file Numerical - Julia\experiments\table6.jl and the resulting table would appear as Numerical - Juliaoutput\Table6.csv
* note that this table contains timing results that can very between dofferent machines and experiments

- Figure 5: this figure requires data both from the Matlab LFA code and the Julia numerical code. 
First run the file LFA - Matlab\experiments\figure5.m and the LFA data will appear as LFA - Matlab\output\Figure5_LFA.m. 
Then run the file Numerical - Julia\experiments\figure5.jl and the resulting data would appear as three files: Numerical - Julia\output\Figure5_Element_numerical.csv, and in the same location Figure5_Plus_numerical.csv and Figure5_RB_numerical.csv. 
Finally run the file Numerical - Julia\experiments\plotFigure5.jl and the plots will be saved as Numerical - Julia\output\Figure5a.png and same for Figure5b.png and Figure5c.png.

- Figure 6: this figure requires data both from the Matlab LFA code and the Julia numerical code. 
First run the file LFA - Matlab\experiments\figure6.m and the LFA data will appear as LFA - Matlab\output\Figure6_LFA.m. 
Then run the file Numerical - Julia\experiments\figure6.jl and the resulting data would appear as three files: Numerical - Julia\output\Figure6_ElementVanka.csv, and in the same location Figure6_PlusVanka.csv and Figure6_RBVanka.csv. 
Finally run the file Numerical - Julia\experiments\plotFigure6.jl and the plots will be saved as Numerical - Julia\output\Figure6a.png and same for Figure6b.png and Figure6c.png.

- Figure 7: run the file Numerical - Julia\experiments\figure7.jl and the resulting data would appear as Numerical - Julia\output\Figure7.csv.
Then run the file Numerical - Julia\experiments\plotFigure6.jl and the plots will be saved as Numerical - Julia\output\Figure7a.png and Figure7b.png.

- Figure 9: run the file Numerical - Julia\experiments\figure9.jl and the resulting data would appear as Numerical - Julia\output\Figure9.csv.
Then run the file Numerical - Julia\experiments\plotFigure9.jl and the plots will be saved as Numerical - Julia\output\Figure9a.png and Figure9b.png.

- Figure 10: run the file Numerical - Julia\experiments\figure10a.jl and then the file  Numerical - Julia\experiments\figure10b.jl and the resulting data would appear as Numerical - Julia\output\Figure10a.csv and in the same location Figure10b_Jacobi.csv and Figure10b_Vanka.csv.
Then run the file Numerical - Julia\experiments\plotFigure10.jl and the plots will be saved as Numerical - Julia\output\Figure10a.png and Figure10b.png.
* note that this figure contains timing results that can very between dofferent machines and experiments

- Figure 12: This experiment require using the overthrust model (SEG), first read the instructions on the end of this document and run the experiment only after the file Numerical - Julia\experiments\BenchmarkModels\3DOverthrust801801187.dat is saved. 
Now, run the file Numerical - Julia\experiments\figure12.jl and the resulting data would appear as Numerical - Julia\output\Figure12.csv.
Then run the file Numerical - Julia\experiments\plotFigure12.jl and the plots will be saved as Numerical - Julia\output\Figure12.png.
* note that this table contains timing results that can very between dofferent machines and experiments



# Overthrust model

Table 6 and Figure 12 require the SEG/EAGE Overthrust velocity model.

The model file is not included in this repository because it is distributed by SEG/EAGE under their own data access terms.

The data can be obtained from:
https://wiki.seg.org/wiki/SEG/EAGE_Salt_and_Overthrust_Models

After downloading the SEG-Y file, convert it using a SEG-Y reader to obtain the velocity values as a one-dimensional array. Reshape the array into a 801×801×187 array and save it as:

Numerical - Julia/experiments/BenchmarkModels/3DOverthrust801801187.dat

The Table 6 and Figure 12 scripts require this file before execution.
