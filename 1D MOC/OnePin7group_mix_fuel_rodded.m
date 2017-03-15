%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [6, 2, 1; % Mixture
    4, 2, 1;  % Fuel
    5, 2, 1]; % CR

input.radii = [ 0.4096, 0.475;
    0.4096, 0.475;
    0.4096, 0.475];
input.pinmesh = [ 15 2 15;
    15 2 15;
    15 2 15];
% Mixtures
input.nmixtures = 1;
input.mixtures = [6, 4, 5];
input.mixvols = [0.5, 0.5];
% Quadrature
input.npol = 8;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 2000;
input.conv_crit = [1.0e-6 1.0e-6];
input.verbose = true;
input.subray = false;

%% Case 1 - 50-50 Mixutre Eigenvalue Case
input.pinmap = [1];