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
input.pinmesh = [ 5 1 5;
    5 1 5;
    5 1 5;
    5 1 5];
% Mixtures
input.nmixtures = 1;
input.mixtures = [6, 1, 5];
input.mixvols = [0.5, 0.5];
% Quadrature
input.npol = 16;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% Solution variables
input.nouters = 2000;
input.conv_crit = [1.0e-5 1.0e-5];
input.verbose = false;
% input.subray = false;
input.subray = true;

%% Base Case - 50-50 Mixture Eignevalue Case
%  mixture, no subray, k-eff = 0.8495354, 251 iterations
%  mixture,    subray, k-eff = 0.9936844, 202 iterations
input.pinmap = [1, 1, 1, 3, 1, 1, 4, 1, 1, 3, 1, 1, 1];
eigenSolver = eigensolverClass(input);
eigenSolver.verbose = true;
eigenSolver.solve();

% %% Case 1 - 50-50 Mixutre FSP
% input.verbose = true;
% fssSolver = FixedSourceSolverClass(input, eigenSolver);
% fssSolver.solve(0, 1);
% names(1) = {sprintf('Mixed')};
% 
% % Case 2 - Control Rod Case
% % pinmap_rodded = 1;
% input.pinmap = [1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
%     1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1, ...
%     1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1];
% fssSolver(2) = FixedSourceSolverClass(input, eigenSolver);
% fssSolver(2).solve(0, 59); %59 Iterations at 1.0e-5
% names(2) = {sprintf('Rodded')};
% 
% % Case 2 - Guide Tube Case
% input.pinmap = [1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
%     1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
%     1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1];
% fssSolver(3) = FixedSourceSolverClass(input, eigenSolver);
% fssSolver(3).solve(0, 110); %110 Iteration as 1.0e-5
% names(3) = {sprintf('Unrodded')};