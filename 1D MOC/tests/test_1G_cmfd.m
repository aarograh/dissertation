function result = test_1G_cmfd( verbose )
%TEST_1G_CMFD Performs a 1-group test with vacuum boundaries and
%with CMFD enabled
%   verbose - Flag to enable/disable output

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pinmap = [1, 2];
input.pitch = 1.0;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [4;
                 6];

input.radii = [ ];
input.pinmesh = [10;
                 10];
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '1group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% Convergence
input.nouters = 1000;
if exist('verbose','var')
    input.verbose = verbose;
else
    input.verbose = true;
end
input.cmfd = true;
%30 - 1.1457170

%% Test Case
solver = eigensolverClass(input);
solver.solve( );

%% Test Solution
ref = 1.2/0.8;
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