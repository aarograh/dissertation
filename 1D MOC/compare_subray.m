function [ ] = compare_subray( input )
%COMPARE_SUBRAY Runs mixed, fully rodded, fully unrodded, and sub-ray
%calculations for the input case.  It is assumed that pin 1 in the pinmap
%is the mixed pin, 2 is the unrodded pin, and 3 is the control rod pin

%% Setup
OnePin7group_mix_fuel_rodded
oldmap = input.pinmap;

%% Eigensolve
eSolver = eigensolverClass(input);
eSolver.solve();

%% Fixed Source Solves
fssMixed = FixedSourceSolverClass(input, eSolver);
fssMixed.solve(0,0);

for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = 3;
    end
end
fssCR = FixedSourceSolverClass(input, eSolver);
fssCR.solve(0,0);

for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = 2;
    end
end
fssNoCR = FixedSourceSolverClass(input, eSolver);
fssNoCR.solve(0,0);

%% Sub-ray fixed source solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = 1;
    end
end
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