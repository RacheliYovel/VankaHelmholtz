# include("../src/SolverFunctions.jl")


# function nnz(A)

#     nnzA = count(!iszero, A)

#     return nnzA

# end


# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                         Table 1                         ")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# n = [64;64;64]

# shift = 0.0
# omega_factor = 1e-6
# order = 4

# levels_arr = [2;3;4;5]

# intergrid_arr = [
#     "tricub",
#     # "mixed",
#     "levdep",
#     "trilin"
# ]

# intergrid_names = [
#     "Tricubic",
#     # "Mixed",
#     "Level-dependent",
#     "Trilinear"
# ]


# results = zeros(
#     Float64,
#     length(intergrid_arr),
#     length(levels_arr)
# )


# for i = 1:length(intergrid_arr)

#     for j = 1:length(levels_arr)

#         levels = levels_arr[j]

#         H, H_s, b,
#         R_arr, P_arr, Ac_arr,
#         LUcoarsest, M_arr, relaxParam =
#             getAcousticHelmholtzMGVankaSetup(
#                 n,
#                 "const",
#                 intergrid_arr[i],
#                 "Vanka",
#                 ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"],
#                 levels,
#                 shift;
#                 order = order,
#                 omega_factor = omega_factor
#             )


#         total_nnz = nnz(H)

#         for k = 1:length(Ac_arr)

#             total_nnz += nnz(Ac_arr[k])

#         end


#         results[i,j] = total_nnz / nnz(H)

#     end

# end


# println()

# println(
#     rpad("Intergrid", 20),
#     " | ",
#     join(
#         [lpad("$levels-level", 10)
#          for levels in levels_arr],
#         " | "
#     )
# )

# println("-"^75)


# for i = 1:length(intergrid_arr)

#     println(
#         rpad(intergrid_names[i], 20),
#         " | ",
#         join(
#             [lpad(string(round(results[i,j], digits=3)), 10)
#              for j = 1:length(levels_arr)],
#             " | "
#         )
#     )

# end


# println()

# println("all done!")



using CSV
using DataFrames

include("../src/SolverFunctions.jl")


function nnz(A)

    nnzA = count(!iszero, A)

    return nnzA

end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 1                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n = [64;64;64]

shift = 0.0
omega_factor = 1e-6
order = 4

levels_arr = [2;3;4;5]

intergrid_arr = [
    "tricub",
    "levdep",
    "trilin"
]

intergrid_names = [
    "Tricubic",
    "Level-dependent",
    "Trilinear"
]


results = zeros(
    Float64,
    length(intergrid_arr),
    length(levels_arr)
)


############################################################
# Compute operator complexities
############################################################

for i = 1:length(intergrid_arr)

    for j = 1:length(levels_arr)

        levels = levels_arr[j]

        H, H_s, b,
        R_arr, P_arr, Ac_arr,
        LUcoarsest, M_arr, relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n,
                "const",
                intergrid_arr[i],
                "Vanka",
                ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"],
                levels,
                shift;
                order = order,
                omega_factor = omega_factor
            )

        total_nnz = nnz(H)

        for k = 1:length(Ac_arr)

            total_nnz += nnz(Ac_arr[k])

        end

        results[i,j] = total_nnz / nnz(H)

    end

end


############################################################
# Save table
############################################################

table = DataFrame(
    Intergrid = intergrid_names,
    Level2 = round.(results[:,1], digits=3),
    Level3 = round.(results[:,2], digits=3),
    Level4 = round.(results[:,3], digits=3),
    Level5 = round.(results[:,4], digits=3)
)

CSV.write(
    "../output/Table1.csv",
    table
)


############################################################
# Display summary
############################################################

println()
println("Table 1 saved to:")
println("../output/Table1.csv")
println()

display(table)

println()
println("all done!")