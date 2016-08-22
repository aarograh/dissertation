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
npol = 1;
% XS Library Info
xsfilename = '1group.xsl';
scattype = 'P0';
% Boundary Conditions
BCond = ['reflecting';'reflecting'];
% BCond = ['vacuum';'vacuum'];
% Convergence
nouters = 100;

%% Test Case
pinmap_rodded = 1;
[solution, mesh] = ...
    MOC_1D(pinmap_rodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond, nouters);

%% Test Solution

ref = 1.2/0.8;
if abs(solution.keff(1) - ref) < 2.0e-6
    display(sprintf('Test Passed! Ref: %g, Test: %g',ref,solution.keff(1)));
else
    display(sprintf('Test Failed! Ref: %g, Test: %g',ref,solution.keff(1)));
end