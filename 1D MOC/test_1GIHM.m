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



% %% Generate Plots
% nfigs = 0;
% 
% % Angular Flux Plots
% ipol=1;
% igroups = 1;
% for j=igroups
%     nfigs = nfigs + 1;
%     figure(nfigs);
%     hold on
%     tmp = 0;
%     for i=1:length(solution)
%         plot(mesh.fsredges,solution(i).angflux(:,ipol,1,j),'linewidth',2)
%         plot(mesh.fsredges,solution(i).angflux(:,ipol,2,j),'linewidth',2)
%         tmp = max(tmp,max(solution(i).angflux(:,ipol,1,j)));
%     end
%     ax = gca;
%     ax.XAxis.TickValues = mesh.xsedges;
%     ax.XAxis.MinorTickValues = mesh.fsredges;
%     ax.XTickLabelRotation = 45;
%     axis([min(mesh.xsedges), max(mesh.xsedges), 0.0, 1.05*tmp]);
%     xlabel('Position (cm)')
%     ylabel('Angular Flux')
%     title('Right-going Angular Flux vs. position')
%     legend('forward','backward')
%     grid on
%     grid minor
% end
% 
% 
% % Scalar Flux Plots
% nfigs = nfigs + 1;
% figure(nfigs);
% hold on;
% cellcenter = 0.5*(mesh.fsredges(1:end-1) + mesh.fsredges(2:end));
% plot(cellcenter,solution(1).scalflux(:,1,1));
% ax = gca;
% ax.XAxis.TickValues = mesh.xsedges;
% ax.XAxis.MinorTickValues = mesh.fsredges;
% axis([min(mesh.xsedges), max(mesh.xsedges), 0.0, 1.05*max(solution.scalflux(:,1,1))]);
% xlabel('Position (cm)')
% ylabel('Scalar Flux (cm^-2)')
% legend('Scalar Flux');
% grid on
% grid minor