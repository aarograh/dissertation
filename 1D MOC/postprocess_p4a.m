
%% Generate Plots
close all;
xgrid = min(fssSolver(1).mesh.xsedges):input.pitch:max(fssSolver(1).mesh.xsedges);
xgridminor = min(fssSolver(1).mesh.xsedges):input.pitch:max(fssSolver(1).mesh.xsedges);
linS = {'-',':','--'};

% Angular Flux Plots
ipol=input.npol/2;
for igroup=1:fssSolver(1).xsLib.ngroups
    figure(igroup);
    hold on;
    tmp = 0;
    tmpflux=zeros(length(fssSolver(1).solution.angflux(1,igroup,ipol,:)),1);
    for i=1:length(fssSolver)
        for j=1:length(fssSolver(i).solution.angflux(1,igroup,ipol,:))
            tmpflux(j)=fssSolver(i).solution.angflux(1,igroup,ipol,j);
        end
        plot(fssSolver(i).mesh.fsredges,tmpflux,'linestyle',linS{i},'linewidth',2);
        tmp = max(tmp,max(fssSolver(i).solution.angflux(1,igroup,ipol,:)));
    end
    ax = gca;
    ax.XAxis.TickValues = xgrid;
    ax.XAxis.MinorTickValues = fssSolver(1).mesh.xsedges;
    ax.XTickLabelRotation = 45;
    ax.FontSize = 11;
    xlabel('Position (cm)');
    ylabel('Angular Flux');
    legend(names);
    grid on;
    grid minor;
    axis([17*input.pitch, 34*input.pitch, -inf, .7*tmp]);
end


% Scalar Flux Plots
figure(8);
hold on;
igroup = 7;
cellcenter = 0.5*(fssSolver(1).mesh.fsredges(1:end-1) + fssSolver(1).mesh.fsredges(2:end));
tmp = 0;
for i=1:length(fssSolver)
    plot(cellcenter,fssSolver(i).solution.scalflux(igroup,:,1),'linestyle',linS{i},'linewidth',2);
    tmp = max(tmp,max(fssSolver(i).solution.scalflux(igroup,:,1)));
end
ax = gca;
ax.XAxis.TickValues = xgrid;
ax.XAxis.MinorTickValues = fssSolver(1).mesh.xsedges;
ax.XTickLabelRotation = 45;
ax.FontSize = 11;
xlabel('Position (cm)')
ylabel('Scalar Flux')
legend(names);
grid on;
grid minor;
axis([17*input.pitch, 34*input.pitch, -inf, .7*tmp]);