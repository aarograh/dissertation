% This script assumes that the inputClass input has been defined
% externally.  It also assumes that an array of pin IDs swappinids
% has been set up.  The first index is the mixed pin, the second
% the absorber pin, and the third the replacement pin
% Also assums that a getref variable has been set to true or false.
% True will cause the fully rodded/unrodded cases to run, while false
% will skip them.  This is used when this script is called repeatedly
% by another script.

%% Setup
oldmap = input.pinmap;

%% Volume-homogenized solve
eSolver(1) = eigensolverClass(input);
eSolver(1).solve();

if getref
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
end

%% Sub-ray fixed source solve
for i=1:length(input.pinmap)
    if oldmap(i) == 1
        input.pinmap(i) = swappinids(1);
    end
end
input.subray = 1;
eSolver(4) = eigensolverClass(input);
eSolver(4).solve();

input.subray = 2;
input.npinSubTrack = 0;
eSolver(5) = eigensolverClass(input);
eSolver(5).solve();
input.npinSubTrack = 1;
eSolver(6) = eigensolverClass(input);
eSolver(6).solve();
input.npinSubTrack = 2;
eSolver(7) = eigensolverClass(input);
eSolver(7).solve();

% Post-process
for i=1:length(eSolver)
    fssSolver(i) = eSolver(i).fss;
end
names = {'Volume-Mixed','Rodded','Unrodded','Subray','Subray - Recombination'};
% postprocess_subrayFSS(input, fssSolver, names);