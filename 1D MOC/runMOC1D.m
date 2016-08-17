close all; clear variables; clc;

%% Rodded Case
% Core Data
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
% Energy Group Info
igroup = 47;
% Solve
[angflux, scalflux, mesh] = ...
    MOC_1D(pinmap, pitch, diag, pinmats, radii, pinmesh, npol, igroup);

%% Generate Plots
hold on

% Angular Flux Plots
figure(1);
ipol=1;
plot(mesh.fsredges,angflux(:,ipol,1),'b','linewidth',2)
ax = gca;
ax.XAxis.TickValues = mesh.xsedges;
ax.XAxis.MinorTickValues = mesh.fsredges;
ax.XTickLabelRotation = 30;
axis([min(mesh.xsedges), max(mesh.xsedges), min(angflux(:,ipol,1)), 1.05*max(angflux(:,ipol,1))]);
xlabel('Position (cm)')
ylabel('Angular Flux')
title('Right-going Angular Flux vs. position')
legend('Rodded')
grid on
grid minor


% Scalar Flux Plots
% cellcenter = 0.5*(mesh.fsredges(1:nfinecells) + mesh.fsredges(2:nfinecells+1));
% figure(2);
% plot(cellcenter,scalflux)
% xlabel('Position (cm)')
% ylabel('Scalar Flux (cm^-2)')