%% Setup
oldmap = input.pinmap;

%% Eigensolve
eSolver = eigensolverClass(input);
eSolver.solve();

%% Fixed Source Solves
% Volume-homogenized solve
fssMixed = FixedSourceSolverClass(input, eSolver);
fssMixed.solve(0,0);

% Fully rodded solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = 3;
    end
end
fssCR = FixedSourceSolverClass(input, eSolver);
fssCR.solve(0,0);

% Fully unrodded solve
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