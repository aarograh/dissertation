function [ ] = postprocess_subrayFSS( input, fssSolver, names)
%POSTPROCESS_SUBRAYFSS Generates plots for fixed source solves to compare

%% Generate Plots
close all;
xgrid = min(fssSolver(1).mesh.xsedges):input.pitch:max(fssSolver(1).mesh.xsedges);
xgridminor = min(fssSolver(1).mesh.xsedges):input.pitch:max(fssSolver(1).mesh.xsedges);
linS = {'-',':','--','-.'};


% Scalar Flux Plots
for igroup=1:fssSolver(1).xsLib.ngroups
	figure(igroup);
	hold on;
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
end

end