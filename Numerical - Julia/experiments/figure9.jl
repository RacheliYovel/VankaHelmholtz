# using PyPlot

# include("../src/SolverFunctions.jl")

# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                         Fig. 9                          ")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# n = [256;256]

# jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
# rb  = ["RB";"RB";"RB";"RB";"RB";"RB"]

# patches_arr = [[jac];[rb]]

# levels = 4

# shifts_arr = [0.17;0.19;0.21;0.23;0.25;0.27]

# models_arr = [
#     "linear025";
#     "wedge025";
#     "linear005";
#     "wedge005"
# ]

# recursive_calls = 1

# ############################################################
# # Compute GMRES iterations
# ############################################################

# results = zeros(
#     Int,
#     length(patches_arr),
#     length(models_arr),
#     length(shifts_arr)
# )

# for i = 1:length(patches_arr)
#     for j = 1:length(models_arr)
#         for k = 1:length(shifts_arr)

#             H,H_s,b,
#             R_arr,P_arr,Ac_arr,
#             LUcoarsest,M_arr,relaxParam =
#                 getAcousticHelmholtzMGVankaSetup(
#                     n,
#                     models_arr[j],
#                     "levdep",
#                     "Vanka",
#                     patches_arr[i],
#                     levels,
#                     shifts_arr[k]
#                 )

#             iterations,elapsed_time =
#                 solveMGVanka(
#                     H,H_s,b,
#                     R_arr,P_arr,Ac_arr,
#                     LUcoarsest,M_arr,
#                     recursive_calls,
#                     levels,
#                     [1;1],
#                     relaxParam,
#                     "Vanka"
#                 )

#             results[i,j,k] = iterations

#         end
#     end
# end

# ############################################################
# # Print tables
# ############################################################

# for i = 1:length(patches_arr)

#     println()

#     if i == 1
#         println("Jacobi")
#     else
#         println("RB")
#     end

#     println()

#     println(
#         rpad("Model",18),
#         " | ",
#         join(
#             [lpad("α = $(shift)",10) for shift in shifts_arr],
#             " | "
#         )
#     )

#     println("-"^105)

#     for j = 1:length(models_arr)

#         println(
#             rpad(models_arr[j],18),
#             " | ",
#             join(
#                 [lpad(string(results[i,j,k]),10)
#                     for k=1:length(shifts_arr)],
#                 " | "
#             )
#         )

#     end

# end

# ############################################################
# # Figure 9a : Jacobi
# ############################################################

# figure()

# for j = 1:length(models_arr)

#     plot(
#         shifts_arr,
#         results[1,j,:],
#         "o-",
#         label=models_arr[j]
#     )

# end

# xlabel("Shift")
# ylabel("Iterations")
# legend()

# savefig("../output/Figure9a.png",
#     dpi=300,
#     bbox_inches="tight")

# ############################################################
# # Figure 9b : RB
# ############################################################

# figure()

# for j = 1:length(models_arr)

#     plot(
#         shifts_arr,
#         results[2,j,:],
#         "o-",
#         label=models_arr[j]
#     )

# end

# xlabel("Shift")
# ylabel("Iterations")
# legend()

# savefig("../output/Figure9b.png",
#     dpi=300,
#     bbox_inches="tight")

# println()
# println("Figure 9 reproduced.")



using CSV
using DataFrames

include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Fig. 9                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n = [256;256]

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
rb  = ["RB";"RB";"RB";"RB";"RB";"RB"]

patches_arr = [[jac];[rb]]

patch_names = [
    "Jacobi",
    "RB"
]

levels = 4

shifts_arr = [0.17;0.19;0.21;0.23;0.25;0.27]

models_arr = [
    "linear025";
    "wedge025";
    "linear005";
    "wedge005"
]

recursive_calls = 1


############################################################
# Compute GMRES iterations
############################################################

results = zeros(
    Int,
    length(patches_arr),
    length(models_arr),
    length(shifts_arr)
)


for i = 1:length(patches_arr)

    for j = 1:length(models_arr)

        for k = 1:length(shifts_arr)

            println(
                patch_names[i],
                ", ",
                models_arr[j],
                ", shift = ",
                shifts_arr[k]
            )

            H,H_s,b,
            R_arr,P_arr,Ac_arr,
            LUcoarsest,M_arr,relaxParam =
                getAcousticHelmholtzMGVankaSetup(
                    n,
                    models_arr[j],
                    "levdep",
                    "Vanka",
                    patches_arr[i],
                    levels,
                    shifts_arr[k]
                )


            iterations,elapsed_time =
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


            results[i,j,k] = iterations

        end
    end
end


############################################################
# Save CSV
############################################################

output = DataFrame(
    smoother = String[],
    model = String[],
    shift = Float64[],
    iterations = Int[]
)


for i = 1:length(patches_arr)
    for j = 1:length(models_arr)
        for k = 1:length(shifts_arr)

            push!(
                output,
                (
                    patch_names[i],
                    models_arr[j],
                    shifts_arr[k],
                    results[i,j,k]
                )
            )

        end
    end
end


CSV.write(
    "../output/Figure9.csv",
    output
)


println()
println("Figure 9 data saved.")
println("all done!")