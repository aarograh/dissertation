close all; clear variables; clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [4, 2, 1; % Fuel
    5, 2, 1; % Control Rod
    1, 2, 1; % Guide Tube
    6, 2, 1]; % 50-50 Mixture

input.radii = [ 0.4096, 0.475;
    0.4, 0.475;
    0.4, 0.475;
    0.4, 0.475];
input.pinmesh = [ 15 2 15;
    15 2 15;
    15 2 15;
    15 2 15];
% Mixtures
input.nmixtures = 1;
input.mixtures = [6, 1, 5];
input.mixvols = [0.5, 0.5];
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 2000;
input.conv_crit = [1.0e-6 1.0e-6];
input.verbose = false;

%% Case 1 - 50-50 Mixutre Eigenvalue Case
input.pinmap = [1, 4, 1, 1, 1, 1, 1, 1, 1, 1];
% solver = MOC_1D(input);
% names(1) = {sprintf('Mixed       - %g',solver(1).solution.keff(1))};
% 
% %% Case 2 - Control Rod Case
% % pinmap_rodded = 1;
% input.pinmap = [1, 2, 1, 1, 1, 1, 1, 1, 1, 1];
% solver(2) = eigensolverClass(input);
% % Copy source, angular flux
% solver(2).solution.fisssrc(:) = solver(1).solution.fisssrc(:);
% solver(2).solution.scalflux(:) = solver(1).solution.scalflux(:);
% solver(2).solution.angflux(:) = solver(1).solution.angflux(:); 
% solver(2).step(true);
% names(2) = {sprintf('Rodded    - %g',solver(1).solution.keff(1))};
% 
% %% Case 3 - Guide Tube Case
% input.pinmap = [1, 3, 1, 1, 1, 1, 1, 1, 1, 1];
% solver(3) = eigensolverClass(input);
% % Copy source, angular flux
% solver(3).solution.fisssrc(:) = solver(1).solution.fisssrc(:);
% solver(3).solution.scalflux(:) = solver(1).solution.scalflux(:);
% solver(3).solution.angflux(:) = solver(1).solution.angflux(:);
% solver(3).step(true);
% names(3) = {sprintf('Unrodded - %g',solver(1).solution.keff(1))};