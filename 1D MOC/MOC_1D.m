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
% Core Information
npins = 3;
pinmap = ['fuel   '; 'control'; 'fuel   '];
pintypes = ['fuel   ';'control'];
cells = [id_mod,id_clad,id_fuel,id_clad,id_mod;
    id_controlmod,id_gtube,id_control,id_gtube,id_controlmod];

% Pin information
pitch = 1.26;
fuelrad = 0.4096;
controlrad = 0.4;
cladrad = [fuelrad, 0.475];
gtuberad = [controlrad, cladrad(2)];
modsize = 1; % 1: narrow water; 2: wide water

% Quadrature
npol = 1;

% Mesh
nmesh_mod = 10;
nmesh_clad = 5;
nmesh_gap = 1;
nmesh_fuel = 30;
nmesh_gtube = nmesh_clad;
nmesh_controlgap = nmesh_gap;
nmesh_controlmod = nmesh_mod;
nmesh_control = 30;

% Cross-sections
ngroups = 47;
igroup = 47;
xstr_list = ones(nmats,ngroups);
xstr_list(id_mod:id_fuel,47) = [4.394922836; 0.299908964; 3.44E-05; 1.489061075];
xstr_list(id_gtube:id_control,47) = [0.299908964; 3.44E-05; 4.39E+00; 19.13440383];
xstr_list(id_controlmod,:) = xstr_list(id_mod,:);
xstr_list(id_controlgap,:) = xstr_list(id_controlgap,:);
source_list = ones(nmats,ngroups);
source_list(id_controlmod,:) = source_list(id_mod,:);
source_list(id_gtube:id_control,47) = [7.58E-04; 1.50E-06; 2.44E-02; 8.44399E-05];
source_list(id_controlgap,:) = source_list(id_gap,:);
source_list(id_mod:id_fuel,47) = [0.136409169; 4.03E-03; 1.50E-06; 6.74E-03];

%% Setup Problem
% Geometry
mesh_list = [nmesh_mod; nmesh_clad; nmesh_gap; nmesh_fuel; nmesh_gtube; nmesh_controlgap; nmesh_controlmod; nmesh_control];
tmp = pitch/2.0;
if (modsize == 2)
    tmp = sqrt(2)*tmp;
end
widths(id_mod,1) = tmp - cladrad(2);
widths(id_fuel,1) = fuelrad*2.0;
widths(id_control,1) = controlrad*2.0;
widths(id_clad,1) = cladrad(2) - cladrad(1);
widths(id_gap,1) = cladrad(1) - fuelrad;
widths(id_gtube,1) = gtuberad(2) - gtuberad(1);
widths(id_controlmod,1) = tmp - gtuberad(2);
widths(id_controlgap,1) = gtuberad(1) - controlrad;
widths = widths./mesh_list;

% Ray Parameters
if (npol == 1)
    mu = pi/4;
    mucos = cos(mu);
    weights = 2*pi;
elseif (npol == 2)
elseif (npol == 4)
elseif (npol == 8)
elseif (npol == 16)
end

% Mesh
[matmesh, finemesh, coarsemesh] = mesh(pintypes, pinmap, cells, mesh_list, widths);
nfinecells = size(finemesh,1)-1;
ncoarsecells = size(coarsemesh,1)-1;

% Setup source and cross-sections
sourcemesh(1:nfinecells,1) = 0.0;
xstrmesh(1:nfinecells,1) = 0.0;
for i=1:nfinecells
    sourcemesh(i) = source_list(matmesh(i),igroup);
    xstrmesh(i) = xstr_list(matmesh(i),igroup);
end

%% Solve Problem
[angflux_edge, angflux_avg, scalflux] = sweep(finemesh, xstrmesh, sourcemesh, 0.0, mucos, weights);

%% Generate Plots
% Setup
hold on
position_center = 0.5*(finemesh(1:nfinecells) + finemesh(2:nfinecells+1));
position_all(1:2:2*nfinecells+1,:,:) = finemesh;
position_all(2:2:2*nfinecells,:,:) = position_center;

% Plot Angular Flux v. position
angflux_plot(1:2:2*nfinecells+1,:,:) = angflux_edge;
angflux_plot(2:2:2*nfinecells,:,:) = angflux_avg;

figure(1);
ipol=1;
plot(position_all,angflux_plot(:,ipol,1),'b','linewidth',2)
plot(position_all,angflux_plot(:,ipol,2),'r','linewidth',2)
ax = gca;
ax.XAxis.TickValues = coarsemesh;
ax.XAxis.MinorTickValues = finemesh;
ax.XTickLabelRotation = 30;
axis([min(coarsemesh), max(coarsemesh), min(angflux_plot(:,ipol,1)), 1.05*max(angflux_plot(:,ipol,1))]);
xlabel('Position (cm)')
ylabel('Angular Flux')
title('Right-going Angular Flux vs. position')
legend('Control Rod','Guide Tube')
grid on
grid minor


% Scalar Flux v. position
% figure(2);
% plot(position_center,scalflux)
% xlabel('Position (cm)')
% ylabel('Scalar Flux (cm^-2)')