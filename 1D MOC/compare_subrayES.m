% This script assumes that the inputClass input has been defined
% externally.  It also assumes that an array of pin IDs swappinids
% has been set up.  The first index is the mixed pin, the second
% the absorber pin, and the third the replacement pin

%% Setup
oldmap = input.pinmap;

%% Volume-homogenized solve
eSolver = eigensolverClass(input);
eSolver.solve();

%% Fully rodded solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(2);
    end
end
eSolver(2) = eigensolverClass(input);
eSolver(2).solve();

%% Fully unrodded solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(3);
    end
end
eSolver(3) = eigensolverClass(input);
eSolver(3).solve();

%% Sub-ray fixed source solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(1);
    end
end
input.subray = true;
eSolver(4) = eigensolverClass(input);
eSolver(4).solve();

% Post-process
fssSolver(1) = eSolver(1).fss;
fssSolver(2) = eSolver(2).fss;
fssSolver(3) = eSolver(3).fss;
fssSolver(4) = eSolver(4).fss;
names = {'Volume-Mixed','Rodded','Unrodded','Subray'};
postprocess_subrayFSS(input, fssSolver, names);