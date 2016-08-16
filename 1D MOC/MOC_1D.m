close all; clear variables; clc;

%% Material IDs
nmats = 8;
id_mod = 1;
id_clad = 2;
id_gap = 3;
id_fuel = 4;
id_gtube = 5;
id_controlgap = 6;
id_controlmod = 7;
id_control = 8;

%% User Options
pinmap = [1, 2, 1];
pitch = 1.26;
diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water

% Pin information
pinmats = [4, 2, 1;
    8, 5, 7];
radii = [ 0.4096, 0.475;
    0.4, 0.475];
pinmesh = [ 20 2 10;
    20 3 10];

% Quadrature
npol = 1;

%% Cross-sections
ngroups = 47;
igroup = 47;
% TODO: Get rid of this, add real cross-sections
xstr_list = genXSData(nmats, ngroups, id_mod, id_clad, id_gap, id_fuel, id_gtube, id_controlgap, ...
    id_controlmod, id_control);
source_list = ones(nmats,ngroups);
source_list(id_controlmod,:) = source_list(id_mod,:);
source_list(id_gtube:id_control,47) = [7.58E-04; 1.50E-06; 2.44E-02; 8.44399E-05];
source_list(id_controlgap,:) = source_list(id_gap,:);
source_list(id_mod:id_fuel,47) = [0.136409169; 4.03E-03; 1.50E-06; 6.74E-03];

%% Quadrature
quad = quadrature(npol);

%% Mesh
[matmesh, finemesh, coarsemesh] = mesh(pinmap, pinmats, radii, pinmesh, pitch, diag);
nfinecells = size(finemesh,1)-1;
ncoarsecells = size(coarsemesh,1)-1;

%% Perform Sweep
[sourcemesh, xstrmesh] = setupFSP(source_list, xstr_list, matmesh, igroup);
[angflux, ~, scalflux] = sweep(finemesh, xstrmesh, sourcemesh, 0.0, quad);

%% Generate Plots
hold on

% Angular Flux Plots
figure(1);
ipol=1;
plot(finemesh,angflux(:,ipol,1),'b','linewidth',2)
ax = gca;
ax.XAxis.TickValues = coarsemesh;
ax.XAxis.MinorTickValues = finemesh;
ax.XTickLabelRotation = 30;
axis([min(coarsemesh), max(coarsemesh), min(angflux(:,ipol,1)), 1.05*max(angflux(:,ipol,1))]);
xlabel('Position (cm)')
ylabel('Angular Flux')
title('Right-going Angular Flux vs. position')
legend('Control Rod')
grid on
grid minor


% Scalar Flux Plots
% cellcenter = 0.5*(finemesh(1:nfinecells) + finemesh(2:nfinecells+1));
% figure(2);
% plot(cellcenter,scalflux)
% xlabel('Position (cm)')
% ylabel('Scalar Flux (cm^-2)')