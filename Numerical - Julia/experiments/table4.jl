include("../src/SolverFunctions.jl")



println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 4                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n_arr = [[[128;128]];[[256;256]];[[512;512]];[[1024;1024]]];
rb = ["RB";"RB";"RB";"RB";"RB";"RB"]
levels_arr = [2;3;4]
shifts_arr = [0.0 0.1 0.25;
              0.0 0.1 0.25;
              0.0 0.1 0.25;
              0.0 0.1 0.3]
models_arr = ["linear025";"wedge025";"linear005";"wedge005"]
recursive_calls = 1

# results[j,k,i] = GMRES iteration count
# j = model
# k = number of levels
# i = grid size

results = zeros(
    Int,
    length(models_arr),
    length(levels_arr),
    length(n_arr)
)

for i = 1:length(n_arr)
    for j = 1:length(models_arr)
        for k = 1:length(levels_arr)

            H, H_s, b,
            R_arr, P_arr, Ac_arr,
            LUcoarsest, M_arr, relaxParam =
                getAcousticHelmholtzMGVankaSetup(
                    n_arr[i],
                    models_arr[j],
                    "levdep",
                    "Vanka",
                    rb,
                    levels_arr[k],
                    shifts_arr[j, k]
                )

            iterations, elapsed_time =
                solveMGVanka(
                    H, H_s, b,
                    R_arr, P_arr, Ac_arr,
                    LUcoarsest, M_arr,
                    recursive_calls,
                    levels_arr[k],
                    [1;1],
                    relaxParam,
                    "Vanka"
                )

            results[j, k, i] = iterations

        end
    end
end



using CSV
using DataFrames


table = DataFrame(
    Grid = String[],
    Frequency = String[],

    Linear025_L2 = String[],
    Linear025_L3 = String[],
    Linear025_L4 = String[],

    Wedge025_L2 = String[],
    Wedge025_L3 = String[],
    Wedge025_L4 = String[],

    Linear005_L2 = String[],
    Linear005_L3 = String[],
    Linear005_L4 = String[],

    Wedge005_L2 = String[],
    Wedge005_L3 = String[],
    Wedge005_L4 = String[]
)



function format_entry(x)

    return string(x)

end



for i = 1:length(n_arr)

    n = n_arr[i][1]

    grid = "$(n)x$(n)"
    freq = "$(n/5)pi"


    push!(
        table,
        (

            grid,
            freq,

            # kappa^2 in [0.25,1]
            format_entry(results[1,1,i]),
            format_entry(results[1,2,i]),
            format_entry(results[1,3,i]),

            format_entry(results[2,1,i]),
            format_entry(results[2,2,i]),
            format_entry(results[2,3,i]),


            # kappa^2 in [0.05,1]
            format_entry(results[3,1,i]),
            format_entry(results[3,2,i]),
            format_entry(results[3,3,i]),

            format_entry(results[4,1,i]),
            format_entry(results[4,2,i]),
            format_entry(results[4,3,i])
        )
    )

end



############################################################
# Save CSV
############################################################

outfile = "../output/Table4.csv"

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