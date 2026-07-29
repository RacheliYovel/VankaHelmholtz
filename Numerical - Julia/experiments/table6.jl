include("../src/SolverFunctions.jl")

using CSV
using DataFrames


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 6                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n_arr = [
    [128;128;40];
    [192;192;60];
    [256;256;80];
    [320;320;100]
]


recursive_calls = 1


############################################################
# Smoother configurations
############################################################

# Jacobi on fine level, Element Vanka afterwards
jac_element = [
    "Jac";
    "Element";
    "Element";
    "Element";
    "Element";
    "Element"
]


# Jacobi on all levels
jacobi = [
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac";
    "Jac"
]


patches_arr = [
    jac_element,
    jacobi
]


patch_names = [
    "Jacobi+Vanka",
    "Jacobi"
]


############################################################
# Storage of raw results
############################################################

raw_table = DataFrame(
    Smoother = String[],
    Grid = String[],
    Levels = Int[],
    Iterations = Int[],
    Time = Float64[]
)



############################################################
# Function to run experiments
############################################################

function run_table6_case!(
    raw_table,
    n_arr_run,
    levels_arr,
    shifts_arr,
    patches_arr,
    patch_names
)

    for i = 1:length(n_arr_run)

        for j = 1:length(levels_arr)

            levels = levels_arr[j]

            for k = 1:length(patches_arr)

                println()
                println(
                    "Running ",
                    patch_names[k],
                    ", grid = ",
                    n_arr_run[i],
                    ", levels = ",
                    levels
                )


                H,H_s,b,
                R_arr,P_arr,Ac_arr,
                LUcoarsest,M_arr,relaxParam =
                    getAcousticHelmholtzMGVankaSetup(
                        n_arr_run[i],
                        "Overthrust",
                        "trilin",
                        "Vanka",
                        patches_arr[k],
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


                grid_string =
                    string(
                        n_arr_run[i][1],
                        "x",
                        n_arr_run[i][2],
                        "x",
                        n_arr_run[i][3]
                    )


                push!(
                    raw_table,
                    (
                        patch_names[k],
                        grid_string,
                        levels,
                        iterations,
                        elapsed_time
                    )
                )

            end

        end

    end

end



############################################################
# Run 3 and 4 level experiments (all grids)
############################################################

run_table6_case!(
    raw_table,
    n_arr,
    [3;4],
    [0.2;0.4],
    patches_arr,
    patch_names
)



############################################################
# Run 5 level experiments (only valid grids)
############################################################

run_table6_case!(
    raw_table,
    [
        [256;256;80];
        [320;320;100]
    ],
    [5],
    [0.5],
    patches_arr,
    patch_names
)



############################################################
# Convert to paper format
############################################################

function format_entry(iter, time)

    return string(
        iter,
        " (",
        round(time, digits=2),
        ")"
    )

end



table = DataFrame(
    Smoother = String[],
    Grid = String[],
    Frequency = String[],
    L3 = String[],
    L4 = String[],
    L5 = String[]
)



frequency_dict = Dict(
    "128x128x40"  => "2.86pi",
    "192x192x60"  => "4.22pi",
    "256x256x80"  => "5.58pi",
    "320x320x100" => "6.97pi"
)



for smoother in patch_names

    for grid in [
        "128x128x40",
        "192x192x60",
        "256x256x80",
        "320x320x100"
    ]

        L3 = "NA"
        L4 = "NA"
        L5 = "NA"


        for row in eachrow(raw_table)

            if row.Smoother == smoother &&
               row.Grid == grid


                if row.Levels == 3

                    L3 = format_entry(
                        row.Iterations,
                        row.Time
                    )


                elseif row.Levels == 4

                    L4 = format_entry(
                        row.Iterations,
                        row.Time
                    )


                elseif row.Levels == 5

                    L5 = format_entry(
                        row.Iterations,
                        row.Time
                    )

                end

            end

        end


        push!(
            table,
            (
                smoother,
                grid,
                frequency_dict[grid],
                L3,
                L4,
                L5
            )
        )

    end


    if smoother != last(patch_names)

        push!(
            table,
            (
                "",
                "",
                "",
                "",
                "",
                ""
            )
        )

    end

end



############################################################
# Save CSV
############################################################

outfile = "../output/Table6.csv"

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