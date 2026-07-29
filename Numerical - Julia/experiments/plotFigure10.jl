using CSV
using DataFrames
using PyPlot


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 10                         ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


############################################################
# Figure 10a
############################################################

fig10a = CSV.read(
    "../output/Figure10a.csv",
    DataFrame
)


figure()

plot(
    fig10a.DOFs,
    fig10a.Jacobi_Trilinear,
    "o-",
    label="Jacobi, Trilinear"
)

plot(
    fig10a.DOFs,
    fig10a.Element_LevelDependent,
    "p-",
    label="Element Vanka, Level-dependent"
)


xlabel("DOFs")
ylabel("Time (s)")
legend()


savefig(
    "../output/Figure10a.png",
    dpi=300,
    bbox_inches="tight"
)



############################################################
# Figure 10b
############################################################

fig10b_jac = CSV.read(
    "../output/Figure10b_Jacobi.csv",
    DataFrame
)

fig10b_vanka = CSV.read(
    "../output/Figure10b_Vanka.csv",
    DataFrame
)


figure()

plot(
    fig10b_jac.Levels_Jacobi,
    fig10b_jac.Time_Jacobi,
    "o-",
    label="Jacobi"
)


plot(
    fig10b_vanka.Levels_Vanka,
    fig10b_vanka.Time_Vanka,
    "p-",
    label="Element Vanka"
)


xlabel("Levels")
ylabel("Time (s)")
legend()


savefig(
    "../output/Figure10b.png",
    dpi=300,
    bbox_inches="tight"
)


println()
println("Figure 10 plots saved.")
println("all done!")