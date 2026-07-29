# using CSV
# using DataFrames

# include("../src/SolverFunctions.jl")


# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                         Fig. 6                          ")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# ############################################################
# # Parameters
# ############################################################

# n = [256;256]

# jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
# element = ["Element";"Element";"Element";"Element";"Element";"Element"]
# plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]
# rb = ["RB";"RB";"RB";"RB";"RB";"RB"]


# patches_arr = [
#     [element];
#     [plus];
#     [rb]
# ]


# smoother_names = [
#     "ElementVanka",
#     "PlusVanka",
#     "RBVanka"
# ]


# levels = 4
# shift = 0.0

# recursive_calls = 1

# tol = 1e-9
# maxit = 100


# ############################################################
# # MG setup and solve
# ############################################################


# for s = 1:length(patches_arr)

#     println()
#     println("Running ", smoother_names[s])


#     H,H_s,b,
#     R_arr,P_arr,Ac_arr,
#     LUcoarsest,M_arr,relaxParam =
#         getAcousticHelmholtzMGVankaSetup(
#             n,
#             "const",
#             "levdep",
#             "Vanka",
#             patches_arr[s],
#             levels,
#             shift
#         )


#     # zero RHS

#     b = 0*b


#     ########################################################
#     # Random complex initial guess
#     ########################################################

#     nodes = size(H,1)

#     guess =
#         rand(nodes) +
#         0.1im*rand(nodes)


#     guess = vec(guess)



#     x,iter,r_arr,x_arr =
#         MGsolver(
#             H_s,
#             b,
#             guess,
#             relaxParam,
#             [1;1],
#             [1;1],
#             levels,
#             recursive_calls,
#             R_arr,
#             P_arr,
#             Ac_arr,
#             LUcoarsest,
#             maxit,
#             tol;
#             relaxType="Vanka",
#             M_arr=M_arr,
#             printres=false
#         )


#     ########################################################
#     # Error history
#     ########################################################

#     err_arr =
#         zeros(length(x_arr))


#     for i=1:length(x_arr)

#         err_arr[i] =
#             norm(x_arr[i]) /
#             norm(x_arr[1])

#     end


#     ########################################################
#     # Save CSV
#     ########################################################

#     output =
#         DataFrame(
#             smoother =
#                 fill(
#                     smoother_names[s],
#                     length(r_arr)
#                 ),
#             iteration =
#                 collect(1:length(r_arr)),
#             residual =
#                 r_arr,
#             error =
#                 err_arr
#         )


#     filename =
#         "../output/Figure6_" *
#         smoother_names[s] *
#         ".csv"


#     CSV.write(
#         filename,
#         output
#     )


#     println("saved ",filename)

# end


# println()
# println("Figure 6 data generated.")


using CSV
using DataFrames
using LinearAlgebra

include("../src/SolverFunctions.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Fig. 6                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Parameters
############################################################

n = [256;256]

element = [
    "Element";"Element";"Element";
    "Element";"Element";"Element"
]

plus = [
    "Plus";"Plus";"Plus";
    "Plus";"Plus";"Plus"
]

rb = [
    "RB";"RB";"RB";
    "RB";"RB";"RB"
]


patches_arr = [
    [element];
    [plus];
    [rb]
]


smoother_names = [
    "ElementVanka",
    "PlusVanka",
    "RBVanka"
]


levels = 2
shift = 0.0

recursive_calls = 1

nu1 = 1
nu2 = 0

tol = 1e-9
maxit = 100


############################################################
# Run MG convergence histories
############################################################

for s = 1:length(patches_arr)

    println()
    println("Running ", smoother_names[s])


    H,H_s,b,
    R_arr,P_arr,Ac_arr,
    LUcoarsest,M_arr,relaxParam =
        getAcousticHelmholtzMGVankaSetup(
            n,
            "const",
            "bicub",
            "Vanka",
            patches_arr[s],
            levels,
            shift
        )


    ########################################################
    # Same setup as old Figure 6:
    # zero RHS + random complex initial guess
    ########################################################

    b = 0*b


    nodes = length(b)

    guess =
        rand(nodes) +
        0.1im*rand(nodes)

    guess = vec(guess)



    ########################################################
    # Multigrid solve
    ########################################################

    x,iter,r_arr,x_arr =
        MGsolver(
            H_s,
            b,
            guess,
            relaxParam,
            nu1,
            nu2,
            levels,
            recursive_calls,
            R_arr,
            P_arr,
            Ac_arr,
            LUcoarsest,
            maxit,
            tol;
            relaxType="Vanka",
            M_arr=M_arr,
            printres=false
        )



    ########################################################
    # Relative error history
    #
    # Same definition as old code:
    # ||x_k|| / ||x_1||
    ########################################################

    err_arr = zeros(length(x_arr))

    for i = 1:length(x_arr)

        err_arr[i] =
            norm(x_arr[i]) /
            norm(x_arr[1])

    end



    ########################################################
    # Save CSV
    ########################################################

    output = DataFrame(
        smoother = fill(
            smoother_names[s],
            length(r_arr)
        ),
        iteration = collect(1:length(r_arr)),
        residual = r_arr,
        error = err_arr
    )


    filename =
        "../output/Figure6_" *
        smoother_names[s] *
        ".csv"


    CSV.write(
        filename,
        output
    )


    println("saved ", filename)

end


println()
println("Figure 6 data generated.")
println("all done!")