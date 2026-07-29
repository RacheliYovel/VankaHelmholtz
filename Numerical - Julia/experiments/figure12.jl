include("../src/SolverFunctions.jl")

using CSV
using DataFrames


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                       Fig. 12                           ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n = [256;256;80]

levels_arr = [3;4;5;6]

shifts_arr = [0.2;0.4;0.5;0.5]

recursive_calls = 1


############################################################
# Three smoother configurations
############################################################

jac = [
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac"
]

element = [
    "Element";
    "Element";
    "Element";
    "Element";
    "Element";
    "Element"
]


# Jacobi on first level, Element Vanka afterwards
jac_element = [
    "Jac";
    "Element";
    "Element";
    "Element";
    "Element";
    "Element"
]


# Keep as a list of patch configurations
# (same order as original Fig. 12 driver)
patches_arr = [
    jac_element,
    jac,
    element
]


patch_names = [
    "Jac+ElementVanka",
    "Jacobi",
    "ElementVanka"
]


############################################################
# Storage
############################################################

table = DataFrame(
    Smoother = String[],
    Levels = Int[],
    Iterations = Int[],
    Time = Float64[]
)


############################################################
# Run
############################################################

for i = 1:length(patches_arr)

    println()
    println("Running ", patch_names[i])

    for j = 1:length(levels_arr)

        levels = levels_arr[j]

        println("  levels = ", levels)


        H,H_s,b,
        R_arr,P_arr,Ac_arr,
        LUcoarsest,M_arr,relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n,
                "Overthrust",
                "trilin",
                "Vanka",
                patches_arr[i],
                levels,
                shifts_arr[j]
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


        push!(
            table,
            (
                patch_names[i],
                levels,
                iterations,
                elapsed_time
            )
        )

    end

end



############################################################
# Save CSV
############################################################

outfile = "../output/Figure12.csv"

CSV.write(
    outfile,
    table
)


println()
println("Saved:")
println(outfile)

display(table)

println()
println("all done!")