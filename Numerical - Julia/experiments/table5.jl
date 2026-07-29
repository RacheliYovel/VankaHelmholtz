# include("../src/SolverFunctions.jl")


# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                         Table 5                         ") 
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# n_arr = [[[48;48;48]];[[64;64;64]];[[96;96;96]];[[128;128;128]]];
# jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
# element = ["Element";"Element";"Element";"Element";"Element";"Element"]
# plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]
# patches_arr = [[jac];[element];[plus]];
# shifts_arr = [0.5;0.4;0.65]
# levels = 5
# recursive_calls = 1

# # results[j,i] = GMRES iteration count
# # j = patch
# # i = grid size

# results = zeros(
#     Int,
#     length(patches_arr),
#     length(n_arr)
# )

# for i = 1:length(n_arr)
#     for j = 1:length(patches_arr)

#         H, H_s, b,
#         R_arr, P_arr, Ac_arr,
#         LUcoarsest, M_arr, relaxParam =
#             getAcousticHelmholtzMGVankaSetup(
#                 n_arr[i],
#                 "const",
#                 "levdep",
#                 "Vanka",
#                 patches_arr[j],
#                 levels,
#                 shifts_arr[j]
#             )

#         iterations, elapsed_time =
#             solveMGVanka(
#                 H, H_s, b,
#                 R_arr, P_arr, Ac_arr,
#                 LUcoarsest, M_arr,
#                 recursive_calls,
#                 levels,
#                 [1;1],
#                 relaxParam,
#                 "Vanka"
#             )

#         results[j, i] = iterations

#     end
# end


# println(
#     rpad("Patch", 18),
#     " | ",
#     join(
#         [lpad("$(n_arr[i][1])³", 10)
#          for i = 1:length(n_arr)],
#         " | "
#     )
# )

# println("-"^75)

# patch_names = [
#     "Jacobi",
#     "Element",
#     "Plus"
# ]

# for j = 1:length(patch_names)

#     println(
#         rpad(patch_names[j], 18),
#         " | ",
#         join(
#             [lpad(string(results[j, i]), 10)
#              for i = 1:length(n_arr)],
#             " | "
#         )
#     )

# end



# println("all done!")




include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 5                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n_arr = [[[48;48;48]];[[64;64;64]];[[96;96;96]];[[128;128;128]]];

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]

patches_arr = [[jac];[element];[plus]];

shifts_arr = [0.5;0.4;0.65]

levels = 5
recursive_calls = 1


############################################################
# results[j,i] = GMRES iteration count
# times[j,i]   = wall-clock time
# j = patch
# i = grid size
############################################################

results = zeros(
    Int,
    length(patches_arr),
    length(n_arr)
)

times = zeros(
    Float64,
    length(patches_arr),
    length(n_arr)
)


############################################################
# Run experiments
############################################################

for i = 1:length(n_arr)

    for j = 1:length(patches_arr)

        H, H_s, b,
        R_arr, P_arr, Ac_arr,
        LUcoarsest, M_arr, relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n_arr[i],
                "const",
                "levdep",
                "Vanka",
                patches_arr[j],
                levels,
                shifts_arr[j]
            )


        iterations, elapsed_time =
            solveMGVanka(
                H, H_s, b,
                R_arr, P_arr, Ac_arr,
                LUcoarsest, M_arr,
                recursive_calls,
                levels,
                [1;1],
                relaxParam,
                "Vanka"
            )


        results[j, i] = iterations
        times[j, i] = elapsed_time

    end

end



############################################################
# Convert to paper-style table
############################################################

using CSV
using DataFrames


function format_entry(j,i)

    return string(
        results[j,i],
        " (",
        round(times[j,i], digits=2),
        ")"
    )

end



table = DataFrame(
    Grid = String[],
    Frequency = String[],
    Jacobi = String[],
    Element = String[],
    Plus = String[]
)



for i = 1:length(n_arr)

    n = n_arr[i][1]

    push!(
        table,
        (
            "$(n)x$(n)x$(n)",
            "$(n/5)*pi",
            format_entry(1,i),
            format_entry(2,i),
            format_entry(3,i)
        )
    )

end



############################################################
# Save CSV
############################################################

outfile = "../output/Table5.csv"

CSV.write(
    outfile,
    table
)


println()
println("Saved:")
println(outfile)

println()

display(table)

println()
println("all done!")