using CSV
using DataFrames

include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Fig. 7                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n = [256;256]

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
rb  = ["RB";"RB";"RB";"RB";"RB";"RB";"RB"]

patches_arr = [[jac];[rb];[rb]]

nu_arr = [[[1;1]];[[1;1]];[[1;0]]]


levels_arr = [2;3;4;5;6;7]


Jacobi11_shift_arr = [0.0;0.1;0.7;2.5;2.5;2.5]
Vanka11_shift_arr  = [0.0;0.1;0.18;0.18;0.18;0.18]
Vanka10_shift_arr  = [0.0;0.08;0.15;0.15;0.15;0.15]


shifts_arr = [
    Jacobi11_shift_arr';
    Vanka11_shift_arr';
    Vanka10_shift_arr'
]


recursive_calls = 1


############################################################
# Compute GMRES iterations
############################################################

results = zeros(
    Int,
    length(patches_arr),
    length(levels_arr)
)


for i = 1:length(levels_arr)

    for j = 1:length(patches_arr)


        println(
            "level = ",
            levels_arr[i],
            ", smoother = ",
            j
        )


        H,H_s,b,
        R_arr,P_arr,Ac_arr,
        LUcoarsest,M_arr,relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n,
                "const",
                "levdep",
                "Vanka",
                patches_arr[j],
                levels_arr[i],
                shifts_arr[j,i]
            )


        iterations,elapsed_time =
            solveMGVanka(
                H,H_s,b,
                R_arr,P_arr,Ac_arr,
                LUcoarsest,M_arr,
                recursive_calls,
                levels_arr[i],
                nu_arr[j],
                relaxParam,
                "Vanka"
            )


        results[j,i] = iterations

    end

end



############################################################
# Save CSV
############################################################

rows = DataFrame(
    smoother = String[],
    level = Int[],
    shift = Float64[],
    iterations = Int[]
)


names = [
    "Jacobi_V11",
    "RB_V11",
    "RB_V10"
]


for j = 1:length(names)

    for i = 1:length(levels_arr)

        push!(
            rows,
            (
                names[j],
                levels_arr[i],
                shifts_arr[j,i],
                results[j,i]
            )
        )

    end

end


CSV.write(
    "../output/Figure7.csv",
    rows
)


println()
println("Figure 7 data saved.")
println("all done!")