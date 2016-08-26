function result = test_3pin7G( verbose )
%TEST1GIHM Performs a 3 pin, 7-group test with a symmetric
% boundary on the right-hand side

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pinmap = [1, 1, 2];
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [4, 3, 2, 1;
                 5, 3, 2, 1];

input.radii = [ 0.4096, 0.42, 0.475;
                0.5, 0.52, 0.6];
input.pinmesh = [10, 1, 1, 3;
                 10, 1, 1, 3];
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 70;
input.verbose = verbose;

%% Test Case
solver = MOC_1D(input);

%% Test Solution
ref = 0.0397620;
if abs(solver.solution.keff(1) - ref) < 1.0e-6 && solver.converged
    display(sprintf('Test Passed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
    result = 1;
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
    result = 0;
end

end