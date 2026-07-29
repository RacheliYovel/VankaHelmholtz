using CSV
using DataFrames
using PyPlot

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                    Plot Fig. 5                          ")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# ----------------------------------------------------------
# Read data
# ----------------------------------------------------------

# Numerical results from Julia
element_num = CSV.read("../output/Figure5_Element_numerical.csv", DataFrame)
plus_num    = CSV.read("../output/Figure5_Plus_numerical.csv", DataFrame)
rb_num      = CSV.read("../output/Figure5_RB_numerical.csv", DataFrame)

# LFA results from MATLAB
lfa = CSV.read("../../LFA - matlab/output/Figure5_LFA.csv", DataFrame)

# ----------------------------------------------------------
# Figure 5a: Element Vanka
# ----------------------------------------------------------

figure()

plot(
    element_num.damping,
    element_num.rho,
    "o-",
    label = "Numerical"
)

plot(
    lfa.damping[lfa.smoother .== "ElementVanka"],
    lfa.rho[lfa.smoother .== "ElementVanka"],
    "s--",
    label = "LFA"
)

xlabel("Damping parameter")
ylabel("Convergence factor")
legend()

savefig("../output/Figure5a.png", dpi=300, bbox_inches="tight")

# ----------------------------------------------------------
# Figure 5b: Plus Vanka
# ----------------------------------------------------------

figure()

plot(
    plus_num.damping,
    plus_num.rho,
    "o-",
    label = "Numerical"
)

plot(
    lfa.damping[lfa.smoother .== "PlusVanka"],
    lfa.rho[lfa.smoother .== "PlusVanka"],
    "s--",
    label = "LFA"
)

xlabel("Damping parameter")
ylabel("Convergence factor")
legend()

savefig("../output/Figure5b.png", dpi=300, bbox_inches="tight")

# ----------------------------------------------------------
# Figure 5c: RB Vanka
# ----------------------------------------------------------

figure()

plot(
    rb_num.damping,
    rb_num.rho,
    "o-",
    label = "Numerical"
)

plot(
    lfa.damping[lfa.smoother .== "RBVanka"],
    lfa.rho[lfa.smoother .== "RBVanka"],
    "s--",
    label = "LFA"
)

xlabel("Damping parameter")
ylabel("Convergence factor")
legend()

savefig("../output/Figure5c.png", dpi=300, bbox_inches="tight")

println()
println("Figure 5 reproduced.")
println("Saved:")
println("  ../output/Figure5a.png")
println("  ../output/Figure5b.png")
println("  ../output/Figure5c.png")