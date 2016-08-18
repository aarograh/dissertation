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
pinmesh = [ 20 3 20;
    20 3 20;
    20 3 20];
% Quadrature
npol = 1;
% XS Library Info
xsfilename = 'c5g7.xsl';
% xsfilename = '2group.xsl';
scattype = 'P0';
% Boundary Conditions
BCond = ['vacuum';'vacuum'];
% Convergence
nouters = 1000;

%% Case 1 - Control Rod Case
% pinmap_rodded = 1;
pinmap_rodded = [1, 2, 1, 1, 1];
[solution(1), mesh] = ...
    MOC_1D(pinmap_rodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond, nouters);
names = {sprintf('Rodded    - %g',solution(1).keff(1))};

%% Case 2 - Guide Tube Case
 pinmap_unrodded = [1, 3, 1, 1, 1];
 [solution(2), ~] = ...
     MOC_1D(pinmap_unrodded, pitch, diag, pinmats, radii, pinmesh, npol, xsfilename, scattype, BCond, nouters);
names(2) = {sprintf('Unrodded - %g',solution(2).keff(1))};

%% Generate Plots

% Angular Flux Plots
ipol=1;
igroups = 1:7;
for j=igroups
    figure(j);
    hold on
    tmp = 0;
    for i=1:length(solution)
        plot(mesh.fsredges,solution(i).angflux(:,ipol,1,j),'linewidth',2)
        tmp = max(tmp,max(solution(i).angflux(:,ipol,1,j)));
    end
    ax = gca;
    ax.XAxis.TickValues = mesh.xsedges;
    ax.XAxis.MinorTickValues = mesh.fsredges;
    ax.XTickLabelRotation = 45;
    axis([min(mesh.xsedges), max(mesh.xsedges), 0.0, 1.05*tmp]);
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
cellcenter = 0.5*(mesh.fsredges(1:end-1) + mesh.fsredges(2:end));
plot(cellcenter,solution(1).scalflux(:,7,1));
plot(cellcenter,solution(2).scalflux(:,7,1));
xlabel('Position (cm)')
ylabel('Scalar Flux (cm^-2)')