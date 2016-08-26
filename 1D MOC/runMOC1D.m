close all; clear variables; clc;

%% General Input Data
% 1: Fuel Pin
% 2: Control Pin
% 3: Guide Tube Pin
input = inputClass();
input.pitch = 1.26;
input.diag = 0; % flat to indicate whether pin moves through narrow (0) or wide (1) water
% Pin information
input.pinmats = [4, 2, 1;
    5, 2, 1;
    1, 2, 1];

input.radii = [ 0.4096, 0.475;
    0.4, 0.475;
    0.4, 0.475];
input.pinmesh = [ 20 3 20;
    20 3 20;
    20 3 20];
% Quadrature
input.npol = 2;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
% xsfilename = '2group.xsl';
input.scattype = 'P0';
% Boundary Conditions
% BCond = ['reflecting';'reflecting'];
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 1000;

%% Case 1 - Control Rod Case
% pinmap_rodded = 1;
input.pinmap = [1, 2, 1, 1, 1, 1, 1, 1, 1, 1];
solver = MOC_1D(input);
names = {sprintf('Rodded    - %g',solver(1).solution.keff(1))};

%% Case 2 - Guide Tube Case
input.pinmap = [1, 3, 1, 1, 1, 1, 1, 1, 1, 1];
solver(2) = MOC_1D(input);
names(2) = {sprintf('Unrodded - %g',solver(2).solution.keff(1))};

%% Generate Plots

% Angular Flux Plots
ipol=1;
igroups = 1:7;
for j=igroups
    figure(j);
    hold on
    tmp = 0;
    for i=1:length(solver)
        plot(solver(i).mesh.fsredges,solver(i).solution.angflux(:,ipol,1,j),'linewidth',2)
        tmp = max(tmp,max(solver(i).solution.angflux(:,ipol,1,j)));
    end
    ax = gca;
    ax.XAxis.TickValues = solver(1).mesh.xsedges;
    ax.XAxis.MinorTickValues = solver(1).mesh.fsredges;
    ax.XTickLabelRotation = 45;
    axis([min(solver(1).mesh.xsedges), max(solver(1).mesh.xsedges), 0.0, 1.05*tmp]);
    xlabel('Position (cm)')
    ylabel('Angular Flux')
    title('Right-going Angular Flux vs. position')
    legend(names)
    grid on
    grid minor
end


% Scalar Flux Plots
figure(8);
hold on;
cellcenter = 0.5*(solver(1).mesh.fsredges(1:end-1) + solver(1).mesh.fsredges(2:end));
plot(cellcenter,solver(1).solution.scalflux(:,7,1));
plot(cellcenter,solver(2).solution.scalflux(:,7,1));
xlabel('Position (cm)')
ylabel('Scalar Flux (cm^-2)')