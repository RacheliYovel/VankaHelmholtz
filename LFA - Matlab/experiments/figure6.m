% Figure 6: LFA convergence factors for optimal damping
%
% This script computes the local convergence factors rho_loc
% corresponding to the optimal real damping parameters used in Figure 6.
%
% The optimization over complex damping parameters was performed separately.
% Since all optimal parameters are real, only the final real values are
% reproduced here.

clear;
close all;
clc;


%% Add the LFA source directory to the MATLAB path

script_dir = fileparts(mfilename("fullpath"));
src_dir = fullfile(script_dir, "..", "src");

addpath(src_dir);


%% Parameters

h = 1/64;
gppw = 10;
gamma = 0.0;

samples = 256;
stencil = 9;
nu = 1;

CGA = "GCA";
intergrid = "cubic";


%% Optimal real damping parameters

smoother = [
    "ElementVanka";
    "PlusVanka";
    "RBVanka"
];

w_opt = [
    0.97;
    0.87;
    0.83
];


%% Compute LFA convergence factors

rho = zeros(length(smoother),1);


for i = 1:length(smoother)

    [rho(i), ~] = AcousticHelmholtzLFA( ...
        h, ...
        gppw, ...
        gamma, ...
        w_opt(i), ...
        smoother(i), ...
        samples, ...
        stencil, ...
        nu, ...
        CGA, ...
        intergrid);

end


%% Save CSV

output_dir = fullfile(script_dir, "..", "output");

if ~exist(output_dir,"dir")
    mkdir(output_dir);
end


output = table( ...
    smoother, ...
    w_opt, ...
    rho, ...
    'VariableNames', ...
    {'smoother','damping','rho'} ...
);


output_file = fullfile(output_dir,"Figure6_LFA.csv");

writetable(output,output_file);


%% Display result

fprintf("\nFigure 6 LFA calculation complete.\n");
fprintf("Results saved to:\n%s\n\n",output_file);

disp(output);