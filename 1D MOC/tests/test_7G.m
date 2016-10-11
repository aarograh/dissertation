function result = test_7G( verbose )
%TEST_7G Performs a 7-group test with vacuum boundaries
%   verbose - Flag to enable/disable output

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
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 95;
if exist('verbose','var')
    input.verbose = verbose;
else
    input.verbose = true;
end

%% Test Case
solver = eigensolverClass(input);
solver.solve( );

%% Test Solution
ref = 0.1384693;
if abs(solver.fss.solution.keff(1) - ref) < 1.0e-6 && solver.converged
    display(sprintf('Test Passed!'));
    result = 1;
elseif abs(solver.fss.solution.keff(1) - ref) < 1.0e-6
    display(sprintf('Result is correct, but did not converge quickly enough!'));
    result = 0;
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.fss.solution.keff(1)));
    result = 0;
end

end