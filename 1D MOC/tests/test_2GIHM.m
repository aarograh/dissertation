function result = test_2GIHM( verbose )
%TEST_2GIHM Performs a 2-group IHM test
%   verbose - Flag to enable/disable output

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pinmap = 1;
input.pitch = 1.0;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = 4;

input.radii = [ ];
input.pinmesh = 1;
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '2group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 60;
if exist('verbose','var')
    input.verbose = verbose;
else
    input.verbose = true;
end

%% Test Case
solver = eigensolverClass(input);
solver.solve( );

%% Setup Reference
clear diag;
xsA = [0.5
0.5];

xsnF = [0.0
1.2];

xsF = [0.0
0.5];

chi = [1.0
0.0];

xsS = [0.5 0.0
0.5 1.0];

xsT = diag(xsA + sum(xsS,1)');

%% Solve
M = xsT - xsS;
phi = M\chi;
ref = xsnF'*phi;

%% Test Solution

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