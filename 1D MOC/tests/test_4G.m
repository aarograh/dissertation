function result = test_4G( verbose )
%TEST1GIHM Performs a 4-group test with vacuum boundaries

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pinmap = 1;
input.pitch = 10.0;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = 1;

input.radii = [ ];
input.pinmesh = 10;
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '4group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 87;
input.verbose = verbose;

%% Test Case
solver = eigensolverClass(input);
solver.solve( );

%% Test Solution
ref = 0.2926313;
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