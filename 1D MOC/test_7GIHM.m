close all; clear variables; %clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
pitch = 1.0;
diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
pinmats = 4;

radii = [ ];
pinmesh = 10;
% Quadrature
npol = 32;
% XS Library Info
xsfilename = 'c5g7.xsl';
scattype = 'P0';
% Boundary Conditions
BCond = ['reflecting';'reflecting'];
% BCond = ['vacuum';'vacuum'];
% Convergence
nouters = 1984;

%% Solve Problem
pinmap_rodded = 1;
solver = ...
    MOC_1D(pinmap_rodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond, nouters);%% Set up XS

%% Setup Reference
clear diag;
xsA = [8.0248E-03
3.7174E-03
2.6769E-02
9.6236E-02
3.0020E-02
1.1126E-01
2.8278E-01];

xsnF = [2.005998E-02
2.027303E-03
1.570599E-02
4.518301E-02
4.334208E-02
2.020901E-01
5.257105E-01];

xsF = [7.21206E-03
8.19301E-04
6.45320E-03
1.85648E-02
1.78084E-02
8.30348E-02
2.16004E-01];

chi = [5.8791E-01
4.1176E-01
3.3906E-04
1.1761E-07
0.0000E+00
0.0000E+00
0.0000E+00];

xsS = [1.27537E-01 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
  4.23780E-02 3.24456E-01 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
  9.43740E-06 1.63140E-03 4.50940E-01 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
  5.51630E-09 3.14270E-09 2.67920E-03 4.52565E-01 1.25250E-04 0.00000E+00 0.00000E+00
  0.00000E+00 0.00000E+00 0.00000E+00 5.56640E-03 2.71401E-01 1.29680E-03 0.00000E+00
  0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.02550E-02 2.65802E-01 8.54580E-03
  0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.00210E-08 1.68090E-02 2.73080E-01];

xsT = diag(xsA + sum(xsS,1)');

%% Solve
M = xsT - xsS;
phi = M\chi;
ref = xsnF'*phi;

%% Test Solution
if abs(solver.solution.keff(1) - ref) < 5.0e-6 && solver.converged
    display(sprintf('Test Passed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
end