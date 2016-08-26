function result = test_4GIHM( verbose )
%TEST1GIHM Performs a 4-group IHM tests
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
input.pinmats = 1;

input.radii = [ ];
input.pinmesh = 1;
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '4group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 855;
input.verbose = verbose;

%% Test Case
solver = ...
    MOC_1D(input);

%% Setup Reference
clear diag;
xsA = [8.0248E-03
3.7174E-03
1.1126E-01
2.8278E-01];

xsnF = [2.005998E-02
2.027303E-03
2.020901E-01
5.257105E-01];

xsF = [7.21206E-03
8.19301E-04
8.30348E-02
2.16004E-01];

chi = [1.0
0.0
0.0000E+00
0.0000E+00];

xsS = [1.27537E-01 0.00000E+00 0.00000E+00 0.00000E+00
4.23780E-02 3.24456E-01 0.00000E+00 0.00000E+00
0.00000E+00 1.00000E-02 2.65802E-01 8.54580E-03
0.00000E+00 1.00000E-03 1.68090E-02 2.73080E-01];

xsT = diag(xsA + sum(xsS,1)');

%% Solve
M = xsT - xsS;
phi = M\chi;
ref = xsnF'*phi;

%% Test Solution
if abs(solver.solution.keff(1) - ref) < 1.0e-6 && solver.converged
    display(sprintf('Test Passed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
    result = 1;
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
    result = 0;
end

end