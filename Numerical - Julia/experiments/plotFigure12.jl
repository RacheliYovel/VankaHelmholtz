using CSV
using DataFrames
using PyPlot


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 12                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Read data
############################################################

data = CSV.read(
    "../output/Figure12.csv",
    DataFrame
)


############################################################
# Plot
############################################################

figure(figsize=(7,5))


smoothers = unique(data.Smoother)


markers = [
    "o",
    "s",
    "^"
]


for i = 1:length(smoothers)

    smoother = smoothers[i]

    subset =
        data[data.Smoother .== smoother, :]


    plot(
        subset.Levels,
        subset.Time,
        marker = markers[i],
        linewidth = 2,
        markersize = 7,
        label = smoother
    )

end


xlabel("number of levels")

ylabel("time (sec)")


xticks(
    sort(unique(data.Levels))
)


grid(true)

legend()


title(
    "Figure 12"
)


############################################################
# Save
############################################################

outfile = "../output/Figure12.png"

savefig(
    outfile,
    dpi = 300,
    bbox_inches = "tight"
)

close()


println()
println("Saved:")
println(outfile)

println()
println("all done!")