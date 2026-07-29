using CSV
using DataFrames

include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Fig. 10b                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n = [128;128;128]

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]

patches_arr = [[jac];[element]]

intergrid_arr = ["levdep";"levdep"]

levels_jac_arr = [3;4;5;6]
levels_vanka_arr = [3;4;5;6;7]

shifts_jac_arr = [0.2;0.5;0.5;0.5]
shifts_vanka_arr = [0.5;0.5;0.5;0.5;0.5]

recursive_calls = 1


############################################################
# Jacobi timings
############################################################

results_jac = zeros(
    Float64,
    length(levels_jac_arr)
)


for i = 1:length(levels_jac_arr)

    levels = levels_jac_arr[i]

    println("Jacobi, levels = ", levels)


    H,H_s,b,
    R_arr,P_arr,Ac_arr,
    LUcoarsest,M_arr,relaxParam =
        getAcousticHelmholtzMGVankaSetup(
            n,
            "linear",
            intergrid_arr[1],
            "Vanka",
            patches_arr[1],
            levels,
            shifts_jac_arr[i]
        )


    iterations, elapsed_time =
        solveMGVanka(
            H,H_s,b,
            R_arr,P_arr,Ac_arr,
            LUcoarsest,M_arr,
            recursive_calls,
            levels,
            [1;1],
            relaxParam,
            "Vanka"
        )


    results_jac[i] = elapsed_time

end


############################################################
# Element Vanka timings
############################################################

results_vanka = zeros(
    Float64,
    length(levels_vanka_arr)
)


for i = 1:length(levels_vanka_arr)

    levels = levels_vanka_arr[i]

    println("Element Vanka, levels = ", levels)


    H,H_s,b,
    R_arr,P_arr,Ac_arr,
    LUcoarsest,M_arr,relaxParam =
        getAcousticHelmholtzMGVankaSetup(
            n,
            "linear",
            intergrid_arr[2],
            "Vanka",
            patches_arr[2],
            levels,
            shifts_vanka_arr[i]
        )


    iterations, elapsed_time =
        solveMGVanka(
            H,H_s,b,
            R_arr,P_arr,Ac_arr,
            LUcoarsest,M_arr,
            recursive_calls,
            levels,
            [1;1],
            relaxParam,
            "Vanka"
        )


    results_vanka[i] = elapsed_time

end


############################################################
# Print table
############################################################

println()

println(
    rpad("Smoother",25),
    " | ",
    join(
        [lpad("L = $L",10) for L in levels_jac_arr],
        " | "
    ),
    " | ",
    lpad("L = 7",10)
)

println("-"^85)


println(
    rpad("Jacobi",25),
    " | ",
    join(
        [lpad(string(round(results_jac[i],digits=3)),10)
         for i=1:length(results_jac)],
        " | "
    )
)


println(
    rpad("Element Vanka",25),
    " | ",
    join(
        [lpad(string(round(results_vanka[i],digits=3)),10)
         for i=1:length(results_vanka)],
        " | "
    )
)


############################################################
# Save CSV
############################################################

output = DataFrame(
    Levels_Jacobi = levels_jac_arr,
    Time_Jacobi = results_jac
)

CSV.write(
    "../output/Figure10b_Jacobi.csv",
    output
)


output = DataFrame(
    Levels_Vanka = levels_vanka_arr,
    Time_Vanka = results_vanka
)

CSV.write(
    "../output/Figure10b_Vanka.csv",
    output
)


println()
println("Figure 10b data saved.")
println("all done!")