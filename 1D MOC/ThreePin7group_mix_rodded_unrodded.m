close all; clear variables; clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Mixtures
input.nmixtures = 1;
input.mixtures = [6, 1, 5];
input.mixvols = [0.5, 0.5];
% Pin information
input.pinmats = [6, 2, 1; % Mixture
    5, 2, 1; % Control Rod
    1, 2, 1; % Guide Tube
    4, 2, 1]; % Fuel

input.radii = [ 0.4, 0.475;
    0.4, 0.475;
    0.4, 0.475;
    0.4096, 0.475];
input.pinmesh = [ 2 1 1;
    5 1 5;
    5 1 5;
    2 1 1];
% Quadrature
input.npol = 2;
% XS Library Info
input.xsfilename = '1group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 2000;
input.conv_crit = [1.0e-5 1.0e-5];
input.verbose = true;

%% Case 1 - 50-50 Mixutre Eigenvalue Case
input.pinmap = [4, 1, 4];

%% set IDs
swappinids = [1,2,3];