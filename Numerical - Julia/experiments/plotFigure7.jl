using CSV
using DataFrames
using PyPlot


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 7                           ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


data = CSV.read(
    "../output/Figure7.csv",
    DataFrame
)


############################################################
# Figure 7a: shift versus levels
############################################################

figure()


for name in unique(data.smoother)

    tmp = data[data.smoother .== name, :]

    plot(
        tmp.level,
        tmp.shift,
        "o-",
        label=name
    )

end


xlabel("Levels")
ylabel("Shift")
legend()


savefig(
    "../output/Figure7a.png",
    dpi=300,
    bbox_inches="tight"
)



############################################################
# Figure 7b: iterations versus levels
############################################################

figure()


for name in unique(data.smoother)

    tmp = data[data.smoother .== name, :]

    plot(
        tmp.level,
        tmp.iterations,
        "o-",
        label=name
    )

end


xlabel("Levels")
ylabel("Iterations")
legend()


savefig(
    "../output/Figure7b.png",
    dpi=300,
    bbox_inches="tight"
)


println("Figure 7 plots saved.")
println("all done!")