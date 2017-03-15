function [ ] = compare_subray( input )
%COMPARE_SUBRAY Runs mixed, fully rodded, fully unrodded, and sub-ray
%calculations for the input case.  It is assumed that pin 1 in the pinmap
%is the mixed pin, 2 is the unrodded pin, and 3 is the control rod pin

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
fssNoCR = FixedSourceSolverClass(input, eSolver);
fssNoCR.solve(0,0);

%% Sub-ray fixed source solve
input.pinmap = 1;
input.subray = true;
fssSubray = FixedSourceSolverClass(input, eSolver);
fssSubray.solve(0,0);

% Post-process
fssSolver(1) = fssMixed;
fssSolver(2) = fssNoCR;
fssSolver(3) = fssCR;
fssSolver(4) = fssSubray;
names = {'Volume-Mixed','Unrodded','Rodded','Subray'};
postprocess_subrayFSS(input, fssSolver, names);

end