function result = test_4GIHM_cmfd( verbose )
%TEST_4GIHM_CMFD Performs a 4-group IHM test with CMFD enabled
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
input.pinmesh = 10;
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = '4group.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% Convergence
input.nouters = 131;
if exist('verbose','var')
    input.verbose = verbose;
else
    input.verbose = true;
end
input.cmfd = true;

%% Test Case
solver = eigensolverClass(input);
solver.solve( );

%% Setup Reference
clear diag;
xsA = [1.17E-02
1.23E-01
1.41E-01
2.83E-01];

xsnF = [2.21E-02
6.09E-02
2.45E-01
5.26E-01];

xsF = [8.03E-03
2.50E-02
1.01E-01
2.16E-01];

chi = [1.00E+00
3.30E-04
0.0000E+00
0.0000E+00];

xsS = [4.94E-01  0.00E+00  0.00E+00  0.00E+00
  1.64E-03  9.06E-01  1.25E-04  0.00E+00
  0.00E+00  5.57E-03  5.49E-01  8.55E-03
  0.00E+00  0.00E+00  1.68E-02  2.73E-01];

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