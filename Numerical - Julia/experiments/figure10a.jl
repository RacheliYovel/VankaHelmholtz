using CSV
using DataFrames

include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Fig. 10a                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n_arr = [[[48;48;48]];[[64;64;64]];[[96;96;96]];[[128;128;128]]]

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]

patches_arr = [[jac];[element]]

intergrid_arr = ["trilin";"levdep"]

shift = 0.5
levels = 5

recursive_calls = 1


############################################################
# Compute elapsed times
#
# results[j,i]
# j = method
# i = grid size
############################################################

results = zeros(
    Float64,
    length(patches_arr),
    length(n_arr)
)


for i = 1:length(n_arr)

    for j = 1:length(patches_arr)

        println(
            "grid = ",
            n_arr[i][1],
            "^3, method = ",
            intergrid_arr[j]
        )

        H,H_s,b,
        R_arr,P_arr,Ac_arr,
        LUcoarsest,M_arr,relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n_arr[i],
                "linear",
                intergrid_arr[j],
                "Vanka",
                patches_arr[j],
                levels,
                shift
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


        results[j,i] = elapsed_time

    end

end


############################################################
# Print table
############################################################

println()

println(
    rpad("Method",25),
    " | ",
    join(
        [
            lpad("$(n_arr[i][1])^3",12)
            for i=1:length(n_arr)
        ],
        " | "
    )
)

println("-"^80)


println(
    rpad("Jacobi Trilinear",25),
    " | ",
    join(
        [
            lpad(string(round(results[1,i],digits=3)),12)
            for i=1:length(n_arr)
        ],
        " | "
    )
)


println(
    rpad("Element Level-dependent",25),
    " | ",
    join(
        [
            lpad(string(round(results[2,i],digits=3)),12)
            for i=1:length(n_arr)
        ],
        " | "
    )
)


############################################################
# Save CSV
#
# x-axis in paper is DOFs = number of grid nodes
############################################################

dofs = [
    (n_arr[i][1]+1)^3
    for i=1:length(n_arr)
]


output = DataFrame(
    DOFs = dofs,
    Jacobi_Trilinear = results[1,:],
    Element_LevelDependent = results[2,:]
)


CSV.write(
    "../output/Figure10a.csv",
    output
)


println()
println("Figure 10a data saved.")
println("all done!")