using CSV
using DataFrames

include("../src/SolverFunctions.jl")
include("../src/MGsolver.jl")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Fig. 5                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# ----------------------------------------------------------
# Parameters
# ----------------------------------------------------------

n = [256;256]

model = "const"
intergrid_type = "bicub"

levels = 2
shift = 0.0

recursive_calls = 1
nu = [1;0]

tol = 1e-9
maxit = 100


damping_arr = Dict(
    "Element" => collect(0.92:0.01:1.00),
    "Plus"    => collect(0.81:0.01:0.89),
    "RB"      => collect(0.79:0.01:0.87)
)


patch_arr = Dict(
    "Element" => ["Element";"Element"],
    "Plus"    => ["Plus";"Plus"],
    "RB"      => ["RB";"RB"]
)



function convergence_factor(r_arr)

    return (r_arr[end]/r_arr[6])^(1/(length(r_arr)-6))

end



for smoother in ["Element","Plus","RB"]

    println()
    println("Smoother: ", smoother)

    damping_values = damping_arr[smoother]

    rho_arr = zeros(length(damping_values))


    for i = 1:length(damping_values)

        w = damping_values[i]

        println("  damping = ", w)


        H,H_s,b,
        R_arr,P_arr,Ac_arr,
        LUcoarsest,M_arr,relaxParam =
            getAcousticHelmholtzMGVankaSetup(
                n,
                model,
                intergrid_type,
                "Vanka",
                patch_arr[smoother],
                levels,
                shift
            )


        # Override the relaxation parameter with the Figure 5 value
        relaxParam[1] = w

        rhs = b
        guess = (0.0 + 0.0*im)*b;


        x, iter, r_arr, x_arr =
            MGsolver(
                H_s,
                rhs,
                guess,
                relaxParam,
                nu[1],
                nu[2],
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


        rho_arr[i] = convergence_factor(r_arr)

    end


    output = DataFrame(
        damping = damping_values,
        rho = rho_arr
    )


    CSV.write(
        "../output/Figure5_" * smoother * "_numerical.csv",
        output
    )

end


println("all done!")