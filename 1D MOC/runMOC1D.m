close all; clear variables; clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
pitch = 1.26;
diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
pinmats = [4, 2, 1;
    5, 2, 1;
    1, 2, 1];

radii = [ 0.4096, 0.475;
    0.4, 0.475;
    0.4, 0.475];
pinmesh = [ 20 2 10;
    20 3 10;
    20 3 10];
% Quadrature
npol = 1;
% XS Library Info
xsfilename = 'c5g7.xsl';
scattype = 'P0';
% Boundary Conditions
BCond = ['vacuum';'vacuum'];

%% Case 1 - Control Rod Case
pinmap_rodded = [1, 2, 1, 1, 1];
[solution(1), mesh] = ...
    MOC_1D(pinmap_rodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond);

%% Case 2 - Guide Tube Case
pinmap_unrodded = [1, 3, 1, 1, 1];
[solution(2), ~] = ...
    MOC_1D(pinmap_unrodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond);

%% Generate Plots
hold on
ipol=1;
igroup = 1;

% Angular Flux Plots
figure(1);
plot(mesh.fsredges,solution(1).angflux(:,ipol,1,igroup),'b','linewidth',2)
plot(mesh.fsredges,solution(2).angflux(:,ipol,1,igroup),'r','linewidth',2)
ax = gca;
ax.XAxis.TickValues = mesh.xsedges;
ax.XAxis.MinorTickValues = mesh.fsredges;
ax.XTickLabelRotation = 45;
tmp = 0;
for i=1:size(solution)
    tmp = max(tmp,solution(i).angflux(:,ipol,1,igroup));
end
axis([min(mesh.xsedges), max(mesh.xsedges), 0.0, 1.05*max(solution(2).angflux(:,ipol,1,igroup))]);
xlabel('Position (cm)')
ylabel('Angular Flux')
title('Right-going Angular Flux vs. position')
legend('Rodded','Unrodded')
grid on
grid minor


% Scalar Flux Plots
% cellcenter = 0.5*(mesh.fsredges(1:nfinecells) + mesh.fsredges(2:nfinecells+1));
% figure(2);
% plot(cellcenter,solution.scalflux)
% xlabel('Position (cm)')
% ylabel('Scalar Flux (cm^-2)')