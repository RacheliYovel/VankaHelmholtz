% Figure 5: LFA convergence factor as a function of the damping parameter
%
% This script computes the LFA convergence factor rho_loc for the three
% Vanka smoothers used in Figure 5:
%
%   (a) Element Vanka: 4-point patch
%   (b) Plus Vanka:    standard 5-point patch
%   (c) RB Vanka:      skew 5-point patch
%
% The resulting data are saved to a CSV file for subsequent combination
% with the numerical convergence factors computed in Julia.

clear;
close all;
clc;

%% Add the LFA source directory to the MATLAB path

script_dir = fileparts(mfilename("fullpath"));
src_dir = fullfile(script_dir, "..", "src");

addpath(src_dir);

%% Parameters used in the LFA calculation

h = 1/64;
gppw = 10;
gamma = 0.0;

samples = 64;
stencil = 9;
nu = 1;
CGA = "GCA";
intergrid = "cubic";

%% Damping parameter ranges used in Figure 5

damping_element = 0.92:0.01:1.00;
damping_plus    = 0.81:0.01:0.89;
damping_rb      = 0.79:0.01:0.87;

%% Calculate the LFA convergence factors

rho_element = zeros(size(damping_element));
rho_plus    = zeros(size(damping_plus));
rho_rb      = zeros(size(damping_rb));

% Figure 5a: Element Vanka

for i = 1:length(damping_element)


w = damping_element(i);

[rho_element(i), ~] = AcousticHelmholtzLFA( ...
    h, ...
    gppw, ...
    gamma, ...
    w, ...
    "ElementVanka", ...
    samples, ...
    stencil, ...
    nu, ...
    CGA, ...
    intergrid);


end

% Figure 5b: Plus Vanka

for i = 1:length(damping_plus)


w = damping_plus(i);

[rho_plus(i), ~] = AcousticHelmholtzLFA( ...
    h, ...
    gppw, ...
    gamma, ...
    w, ...
    "PlusVanka", ...
    samples, ...
    stencil, ...
    nu, ...
    CGA, ...
    intergrid);


end

% Figure 5c: RB Vanka

for i = 1:length(damping_rb)


w = damping_rb(i);

[rho_rb(i), ~] = AcousticHelmholtzLFA( ...
    h, ...
    gppw, ...
    gamma, ...
    w, ...
    "RBVanka", ...
    samples, ...
    stencil, ...
    nu, ...
    CGA, ...
    intergrid);


end

%% Save the results in a single CSV file

damping = [
damping_element(:);
damping_plus(:);
damping_rb(:)
];

rho = [
rho_element(:);
rho_plus(:);
rho_rb(:)
];

smoother = [
repmat("ElementVanka", length(damping_element), 1);
repmat("PlusVanka",    length(damping_plus),    1);
repmat("RBVanka",      length(damping_rb),      1)
];

figure_panel = [
repmat("5a", length(damping_element), 1);
repmat("5b", length(damping_plus),    1);
repmat("5c", length(damping_rb),      1)
];

output = table( ...
figure_panel, ...
smoother, ...
damping, ...
rho, ...
'VariableNames', {'figure_panel', 'smoother', 'damping', 'rho'});
output_dir = fullfile(script_dir, "..", "output");

if ~exist(output_dir, "dir")
mkdir(output_dir);
end

output_file = fullfile(output_dir, "Figure5_LFA.csv");

writetable(output, output_file);

%% Display a summary

fprintf("\n");
fprintf("Figure 5 LFA calculation complete.\n");
fprintf("Results saved to:\n");
fprintf("%s\n\n", output_file);

disp(output);
