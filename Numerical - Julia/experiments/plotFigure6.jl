# using CSV
# using DataFrames
# using PyPlot


# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                    Plot Fig. 6                          ")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# ############################################################
# # Read MATLAB LFA data
# ############################################################

# lfa =
#     CSV.read(
#         "../output/Figure6_LFA.csv",
#         DataFrame
#     )


# ############################################################
# # Plot one smoother
# ############################################################

# function plot_panel(
#     smoother,
#     filename,
#     title
# )

#     julia_data =
#         CSV.read(
#             "../output/Figure6_" *
#             smoother *
#             ".csv",
#             DataFrame
#         )


#     rho =
#         lfa.rho[
#             lfa.smoother .== smoother
#         ][1]


#     t =
#         julia_data.iteration


#     figure()


#     semilogy(
#         t,
#         julia_data.residual,
#         label="Relative residual"
#     )


#     semilogy(
#         t,
#         julia_data.error,
#         label="Relative error"
#     )


#     semilogy(
#         t,
#         rho .^ t,
#         "--",
#         label="Theoretical prediction"
#     )


#     xlabel("iterations")

#     ylabel(
#         L"\|r^{(k)}\|/\|r^{(0)}\|"
#     )


#     ylim(
#         1e-10,
#         1
#     )


#     legend()

#     title(title)


#     savefig(
#         "../output/" * filename,
#         dpi=300,
#         bbox_inches="tight"
#     )

#     close()

# end



# ############################################################
# # Generate panels
# ############################################################


# plot_panel(
#     "ElementVanka",
#     "Figure6a.png",
#     "Element Vanka"
# )


# plot_panel(
#     "PlusVanka",
#     "Figure6b.png",
#     "Plus Vanka"
# )


# plot_panel(
#     "RBVanka",
#     "Figure6c.png",
#     "RB Vanka"
# )


# println()
# println("Figure 6 figures saved.")
# println("all done!")



using CSV
using DataFrames
using PyPlot


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 6                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Directories
############################################################

script_dir = @__DIR__

julia_output_dir = normpath(
    joinpath(script_dir, "..", "output")
)

matlab_lfa_file = normpath(
    joinpath(script_dir, "..", "..", "LFA - matlab", "output", "Figure6_LFA.csv")
)


############################################################
# Read MATLAB LFA data
############################################################

if !isfile(matlab_lfa_file)

    error(
        "Missing MATLAB LFA file:\n",
        matlab_lfa_file
    )

end


lfa = CSV.read(
    matlab_lfa_file,
    DataFrame
)



############################################################
# Plot one smoother
############################################################

function plot_panel(
    smoother,
    filename,
    plot_title,
    julia_output_dir,
    lfa
)


    julia_file = joinpath(
        julia_output_dir,
        "Figure6_" * smoother * ".csv"
    )


    if !isfile(julia_file)

        error(
            "Missing Julia data file:\n",
            julia_file
        )

    end


    julia_data = CSV.read(
        julia_file,
        DataFrame
    )


    rho =
        lfa.rho[
            lfa.smoother .== smoother
        ][1]


    t = julia_data.iteration


    figure()


    semilogy(
        t,
        julia_data.residual,
        label="Relative residual"
    )


    semilogy(
        t,
        julia_data.error,
        label="Relative error"
    )


    semilogy(
        t,
        rho .^ t,
        "--",
        label="Theoretical prediction"
    )


    xlabel("iterations")

    ylabel(
        L"\|r^{(k)}\|/\|r^{(0)}\|"
    )


    ylim(
        1e-10,
        1
    )


    legend()

    title(plot_title)


    savefig(
        joinpath(
            julia_output_dir,
            filename
        ),
        dpi=300,
        bbox_inches="tight"
    )


    close()

end



############################################################
# Generate Figure 6 panels
############################################################

plot_panel(
    "ElementVanka",
    "Figure6a.png",
    "Element Vanka",
    julia_output_dir,
    lfa
)


plot_panel(
    "PlusVanka",
    "Figure6b.png",
    "Plus Vanka",
    julia_output_dir,
    lfa
)


plot_panel(
    "RBVanka",
    "Figure6c.png",
    "RB Vanka",
    julia_output_dir,
    lfa
)


println()
println("Figure 6 figures saved.")
println("all done!")