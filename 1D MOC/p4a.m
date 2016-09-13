close all; clear variables; clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [4, 2, 1; % Fuel
    5, 2, 1; % Control Rod
    1, 2, 1; % Guide Tube
    6, 2, 1]; % 50-50 Mixture

input.radii = [ 0.4096, 0.475;
    0.4, 0.475;
    0.4, 0.475;
    0.4, 0.475];
input.pinmesh = [ 15 2 15;
    15 2 15;
    15 2 15;
    15 2 15];
% Quadrature
input.npol = 16;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['reflecting';'reflecting'];
% Convergence
input.nouters = 2000;
input.conv_crit = [1.0e-5 1.0e-5];

%% Case 1 - 50-50 Mixutre Eigenvalue Case
input.pinmap = [1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
    1, 1, 4, 1, 1, 4, 1, 1, 3, 1, 1, 4, 1, 1, 4, 1, 1, ...
    1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1];
solver = MOC_1D(input);
names(1) = {sprintf('Mixed       - %g',solver(1).solution.keff(1))};

%% Case 2 - Control Rod Case
% pinmap_rodded = 1;
input.pinmap = [1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
    1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1, ...
    1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1];
solver(2) = eigensolverClass(input);
% Copy source, angular flux
solver(2).solution.fisssrc(:) = solver(1).solution.fisssrc(:);
solver(2).solution.scalflux(:) = solver(1).solution.scalflux(:);
solver(2).solution.angflux(:) = solver(1).solution.angflux(:); 
solver(2).step(true);
names(2) = {sprintf('Rodded    - %g',solver(1).solution.keff(1))};

%% Case 2 - Guide Tube Case
input.pinmap = [1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
    1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, ...
    1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1, 3, 1, 1];
solver(3) = eigensolverClass(input);
% Copy source, angular flux
solver(3).solution.fisssrc(:) = solver(1).solution.fisssrc(:);
solver(3).solution.scalflux(:) = solver(1).solution.scalflux(:);
solver(3).solution.angflux(:) = solver(1).solution.angflux(:);
solver(3).step(true);
names(3) = {sprintf('Unrodded - %g',solver(1).solution.keff(1))};

%% Generate Plots

% Angular Flux Plots
ipol=input.npol/2;
for igroup=1:solver(1).xsLib.ngroups
    figure(igroup);
    hold on;
    tmp = 0;
    for i=1:length(solver)
        plot(solver(i).mesh.fsredges,solver(i).solution.angflux(:,ipol,1,igroup),'linewidth',2);
        tmp = max(tmp,max(solver(i).solution.angflux(:,ipol,1,igroup)));
    end
    ax = gca;
    ax.XAxis.TickValues = solver(1).mesh.xsedges;
    ax.XAxis.MinorTickValues = solver(1).mesh.fsredges;
    ax.XTickLabelRotation = 45;
    axis([min(solver(1).mesh.xsedges), max(solver(1).mesh.xsedges), 0.0, 1.05*tmp]);
    xlabel('Position (cm)');
    ylabel('Angular Flux');
    title(sprintf('Right-going Angular Flux vs. Position, Group %i',igroup));
    legend(names);
    grid on;
end


% Scalar Flux Plots
figure(8);
hold on;
igroup = 7;
cellcenter = 0.5*(solver(1).mesh.fsredges(1:end-1) + solver(1).mesh.fsredges(2:end));
tmp = 0;
for i=1:length(solver)
    plot(cellcenter,solver(i).solution.scalflux(:,igroup,1),'linewidth',2);
    tmp = max(tmp,max(solver(i).solution.scalflux(:,igroup,1)));
end
ax = gca;
ax.XAxis.TickValues = solver(1).mesh.xsedges;
ax.XAxis.MinorTickValues = solver(1).mesh.fsredges;
ax.XTickLabelRotation = 45;
axis([min(solver(1).mesh.xsedges), max(solver(1).mesh.xsedges), 0.0, 1.05*tmp]);
xlabel('Position (cm)')
ylabel('Scalar Flux')
title(sprintf('Scalar Flux vs. Position, Group %i',igroup));
legend(names);
grid on;