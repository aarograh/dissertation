% This script assumes that the inputClass input has been defined
% externally.  It also assumes that an array of pin IDs swappinids
% has been set up.  The first index is the mixed pin, the second
% the absorber pin, and the third the replacement pin

%% Setup
oldmap = input.pinmap;
ninners = 0;

%% Eigensolve
eSolver = eigensolverClass(input);
eSolver.solve();

%% Fixed Source Solves
% Volume-homogenized solve
fssMixed = FixedSourceSolverClass(input, eSolver);
fssMixed.solve(0,ninners);

% Fully rodded solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(2);
    end
end
fssCR = FixedSourceSolverClass(input, eSolver);
fssCR.solve(0,ninners);

% Fully unrodded solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(3);
    end
end
fssNoCR = FixedSourceSolverClass(input, eSolver);
fssNoCR.solve(0,ninners);

%% Sub-ray fixed source solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(1);
    end
end
input.subray = true;
fssSubray = FixedSourceSolverClass(input, eSolver);
fssSubray.solve(0,ninners);

% Post-process
fssSolver(1) = fssMixed;
fssSolver(2) = fssNoCR;
fssSolver(3) = fssCR;
fssSolver(4) = fssSubray;
names = {'Volume-Mixed','Rodded','Unrodded','Subray'};
postprocess_subrayFSS(input, fssSolver, names);