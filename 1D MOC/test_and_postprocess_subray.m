close all;
clear all;
clc;

%% Setup
OnePin7group_mix_fuel_rodded

%% Eigensolve
eSolver = eigensolverClass(input);
eSolver.solve();

%% Fixed Source Solves
fssMixed = FixedSourceSolverClass(input, eSolver);
fssMixed.solve(0,0);

input.pinmap = 3;
fssCR = FixedSourceSolverClass(input, eSolver);
fssCR.solve(0,0);

input.pinmap = 2;
fssFuel = FixedSourceSolverClass(input, eSolver);
fssFuel.solve(0,0);

%% Sub-ray fixed source solve
input.pinmap = 1;
input.subray = true;
fssSubray = FixedSourceSolverClass(input, eSolver);
fssSubray.solve(0,0);

% Post-process
fssSolver(1) = fssMixed;
fssSolver(2) = fssFuel;
fssSolver(3) = fssCR;
fssSolver(4) = fssSubray;
names = {'Volume-Homogenized','Fuel','Control Rod','Sub-Ray'};
postprocess_subray