close all; clear variables; %clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
pitch = 10.0;
diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
pinmats = 1;

radii = [ ];
pinmesh = 10;
% Quadrature
npol = 32;
% XS Library Info
xsfilename = '4group.xsl';
scattype = 'P0';
% Boundary Conditions
BCond = ['vacuum';'vacuum'];
% Convergence
nouters = 89;

%% Test Case
pinmap_rodded = 1;
solver = ...
    MOC_1D(pinmap_rodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond, nouters);

%% Test Solution
ref = 0.2926313;
if abs(solver.solution.keff(1) - ref) < 1.0e-6 && solver.converged
    display(sprintf('Test Passed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
else
    display(sprintf('Test Failed! Ref: %0.7f, Test: %0.7f',ref,solver.solution.keff(1)));
end