include("../src/SolverFunctions.jl")

using CSV
using DataFrames


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 3                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n_arr = [[[128;128]];[[256;256]];[[512;512]]]

jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]
rb = ["RB";"RB";"RB";"RB";"RB";"RB"]

patches_arr = [
    [jac];
    [element];
    [plus];
    [rb]
]

patch_names = [
    "Jacobi, alpha=0.3",
    "Element, alpha=0.25",
    "Plus, alpha=0.25",
    "RB, alpha=0.18"
]


shifts_arr = [0.3;0.25;0.25;0.18]

intergrid_arr = [
    "bicub",
    "mixed",
    "levdep"
]

inners = [5;100]

levels = 4
recursive_calls = 2



# results[inner, patch, grid, intergrid]

results = zeros(
    Int,
    length(inners),
    length(patches_arr),
    length(n_arr),
    length(intergrid_arr)
)



############################################################
# Original Table 3 computation
############################################################

for i = 1:length(n_arr)

    for j = 1:length(patches_arr)

        for k = 1:length(intergrid_arr)

            for l = 1:length(inners)


                H, H_s, b,
                R_arr, P_arr, Ac_arr,
                LUcoarsest, M_arr, relaxParam =
                    getAcousticHelmholtzMGVankaSetup(
                        n_arr[i],
                        "const",
                        intergrid_arr[k],
                        "Vanka",
                        patches_arr[j],
                        levels,
                        shifts_arr[j]
                    )


                iter =
                    solveMGVanka(
                        H, H_s, b,
                        R_arr, P_arr, Ac_arr,
                        LUcoarsest, M_arr,
                        recursive_calls,
                        levels,
                        [1;1],
                        relaxParam,
                        "Vanka";
                        inner = inners[l]
                    )


                # current solveMGVanka returns (iterations, ...)
                results[l,j,i,k] = iter[1]


            end
        end
    end
end


############################################################
# Convert to paper-style table
############################################################

function format_entry(j,i,k)

    return string(
        results[1,j,i,k],
        " (",
        results[2,j,i,k],
        ")"
    )

end


table = DataFrame(
    Smoother = String[],
    Grid = String[],
    Frequency = String[],
    Bicubic = String[],
    Mixed = String[],
    LevDep = String[]
)


for j = 1:length(patches_arr)

    for i = 1:length(n_arr)

        n = n_arr[i][1]

        push!(
            table,
            (
                patch_names[j],
                "$(n)x$(n)",
                "$(n/5)*pi",
                format_entry(j,i,1),
                format_entry(j,i,2),
                format_entry(j,i,3)
            )
        )

    end


    if j == 2

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

outfile = "../output/Table3.csv"

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