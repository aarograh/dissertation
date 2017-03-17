close all;
clear all;
clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Mixtures
input.nmixtures = 1;
input.mixtures = [6, 4, 5];
input.mixvols = [0.5, 0.5];
% Pin information
input.pinmats = [6, 2, 1; % Mixture
    5, 2, 1;  % CR
    4, 2, 1]; % Fuel

input.radii = [ 0.4096, 0.475;
    0.4096, 0.475;
    0.4096, 0.475];
input.pinmesh = [ 5 1 5;
    5 1 5;
    5 1 5];
% Quadrature
input.npol = 16;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% Convergence
input.nouters = 500;
input.conv_crit = [1.0e-5 1.0e-5];
input.verbose = true;

%% Case 1 - 50-50 Mixutre Eigenvalue Case
input.pinmap = [1];

%% set IDs
swappinids = [1,2,3];