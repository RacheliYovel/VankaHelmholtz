

# using CSV
# using DataFrames
# using PyPlot


# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("                    Plot Fig. 9                          ")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# data = CSV.read(
#     "../output/Figure9.csv",
#     DataFrame
# )


# models = [
#     "linear025",
#     "wedge025",
#     "linear005",
#     "wedge005"
# ]


# ############################################################
# # Plot settings
# ############################################################

# # The paper shows:
# # Figure 9a = RB
# # Figure 9b = Jacobi
# #
# # Very large iteration counts are outside the visible range,
# # so they are not displayed.

# max_iterations = 100


# ############################################################
# # Function for plotting one smoother
# ############################################################

# function make_plot(smoother, filename, plot_title)

#     figure()

#     for model in models

#         mask = (
#             (data.smoother .== smoother) .&
#             (data.model .== model) .&
#             (data.iterations .<= max_iterations)
#         )

#         sub = data[mask, :]


#         plot(
#             sub.shift,
#             sub.iterations,
#             "o-",
#             label=model
#         )

#     end


#     xlabel("Shift")
#     ylabel("Iterations")

#     ylim(0, max_iterations)

#     title(plot_title)

#     legend()


#     savefig(
#         "../output/" * filename,
#         dpi=300,
#         bbox_inches="tight"
#     )

#     close()

# end


# ############################################################
# # Generate Figure 9 panels
# ############################################################

# # Figure 9a in the paper: RB smoother
# make_plot(
#     "RB",
#     "Figure9a.png",
#     "RB smoother"
# )


# # Figure 9b in the paper: Jacobi smoother
# make_plot(
#     "Jacobi",
#     "Figure9b.png",
#     "Jacobi smoother"
# )


# println()
# println("Figures 9a and 9b saved.")
# println("all done!")


using CSV
using DataFrames
using PyPlot


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 9                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


data = CSV.read(
    "../output/Figure9.csv",
    DataFrame
)


models = [
    "linear025",
    "wedge025",
    "linear005",
    "wedge005"
]


############################################################
# Function for plotting one smoother
############################################################

function make_plot(smoother, filename, plot_title, max_iterations)

    figure()

    for model in models

        mask = (
            (data.smoother .== smoother) .&
            (data.model .== model) .&
            (data.iterations .<= max_iterations)
        )

        sub = data[mask, :]


        plot(
            sub.shift,
            sub.iterations,
            "o-",
            label=model
        )

    end


    xlabel("Shift")
    ylabel("Iterations")

    ylim(0, max_iterations)

    title(plot_title)

    legend()


    savefig(
        "../output/" * filename,
        dpi=300,
        bbox_inches="tight"
    )

    close()

end


############################################################
# Generate Figure 9 panels
############################################################

# Figure 9a in the paper: RB smoother
make_plot(
    "RB",
    "Figure9a.png",
    "RB smoother",
    250
)


# Figure 9b in the paper: Jacobi smoother
make_plot(
    "Jacobi",
    "Figure9b.png",
    "Jacobi smoother",
    500
)


println()
println("Figures 9a and 9b saved.")
println("all done!")