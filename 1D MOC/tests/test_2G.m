function result = test_2G( verbose )
%TEST1GIHM Performs a 2-group test with vacuum boundaries

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pinmap = 1;
input.pitch = 10.0;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = 4;

input.radii = [ ];
input.pinmesh = 10;
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '2group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 119;
input.verbose = verbose;

%% Test Case
solver = MOC_1D(input);

%% Test Solution
ref = 1.1185010;
if abs(solver.solution.keff(1) - ref) < 1.0e-6 && solver.converged
    display(sprintf('Test Passed!'));
    result = 1;
elseif abs(solver.solution.keff(1) - ref)
    display(sprintf('Result is correct, but did not converge quickly enough!'));
    result = 0;
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
    result = 0;
end

end