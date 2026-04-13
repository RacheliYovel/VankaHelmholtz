h = 1/64;
gppw = 10; % grid points per wavelength
gamma = 0.0;
shift = 0.0; % complex shift of CSLP
omega = (2 * pi)/(h * gppw);
gamma = gamma + shift*omega;
stencil = 9; % 5 point 2nd order stencil of 9 point 4th order stencils
% samples = 256;
samples = 64;
nu = 1;
CGA = "GCA";
% intergrid = "FW";
% intergrid = "mixed";
intergrid = "cubic";

smoother_arr = ["Jacobi";"ElementVanka";"PlusVanka";"RBVanka"];

for i=1:length(smoother_arr)
    [rho_loc,mu_loc,w] = optimizeDamping(h,gppw,gamma,smoother_arr(i),samples,stencil,nu,CGA,intergrid);
    disp("-------------------------")
    disp("smoother:")
    disp(smoother_arr(i))
    disp("rho_loc:")
    disp(rho_loc)
    disp("mu_loc:")
    disp(mu_loc)
    disp("damping parameter:")
    disp(w)
    disp("-------------------------")
end