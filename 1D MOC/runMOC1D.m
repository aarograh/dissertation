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
input.pinmesh = [ 20 3 20;
    20 3 20;
    20 3 20;
    20 3 20];
% Quadrature
input.npol = 32;
% XS Library Info
input.xsfilename = 'c5g7.xsl';
input.scattype = 'P0';
% Boundary Conditions
input.BCond = ['vacuum';'vacuum'];
% Convergence
input.nouters = 200;
input.conv_crit = [1.0e-6 1.0e-6];

%% Base Case - 50-50 Mixutre Eigenvalue Case
input.pinmap = [1, 4, 1];
eigenSolver = eigensolverClass(input);
eigenSolver.verbose = true;
eigenSolver.solve();

%% Case 1 - 50-50 FSP
input.verbose = true;
fssSolver = FixedSourceSolverClass(input, eigenSolver);
fssSolver.solve(0, 1); %110 Iterations at 1.0e-5
names(1) = {sprintf('Mixed')};

%% Case 2 - Control Rod FSP
% pinmap_rodded = 1;
input.pinmap = [1, 2, 1];
fssSolver(2) = FixedSourceSolverClass(input, eigenSolver);
fssSolver(2).solve(0, 1); %59 Iterations at 1.0e-5
names(2) = {sprintf('Rodded')};

%% Case 2 - Guide Tube FSP
input.pinmap = [1, 3, 1];
fssSolver(3) = FixedSourceSolverClass(input, eigenSolver);
fssSolver(3).solve(0, 1); %110 Iteration as 1.0e-5
names(3) = {sprintf('Unrodded')};

%% Generate Plots

% Angular Flux Plots
ipol=1;
for igroup=1:fssSolver(1).xsLib.ngroups
    figure(igroup);
    hold on;
    tmp = 0;
    for i=1:length(fssSolver)
        if i == 1
            width = 4;
        else
            width = 2;
        end
        plot(fssSolver(i).mesh.fsredges,fssSolver(i).solution.angflux(:,ipol,1,igroup),'linewidth',width);
        tmp = max(tmp,max(fssSolver(i).solution.angflux(:,ipol,1,igroup)));
    end
    ax = gca;
    ax.XAxis.TickValues = fssSolver(1).mesh.xsedges;
    ax.XAxis.MinorTickValues = fssSolver(1).mesh.fsredges;
    ax.XTickLabelRotation = 45;
    axis([min(fssSolver(1).mesh.xsedges), max(fssSolver(1).mesh.xsedges), 0.0, 1.05*tmp]);
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
cellcenter = 0.5*(fssSolver(1).mesh.fsredges(1:end-1) + fssSolver(1).mesh.fsredges(2:end));
tmp = 0;
for i=1:length(fssSolver)
        if i == 1
            width = 4;
        else
            width = 2;
        end
    plot(cellcenter,fssSolver(i).solution.scalflux(:,igroup,1),'linewidth',width);
    tmp = max(tmp,max(fssSolver(i).solution.scalflux(:,igroup,1)));
end
ax = gca;
ax.XAxis.TickValues = fssSolver(1).mesh.xsedges;
ax.XAxis.MinorTickValues = fssSolver(1).mesh.fsredges;
ax.XTickLabelRotation = 45;
axis([min(fssSolver(1).mesh.xsedges), max(fssSolver(1).mesh.xsedges), 0.0, 1.05*tmp]);
xlabel('Position (cm)')
ylabel('Scalar Flux')
title(sprintf('Scalar Flux vs. Position, Group %i',igroup));
legend(names);
grid on;