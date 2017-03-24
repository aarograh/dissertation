close all;
linS = {'-',':','--','-.'};
fssIndexes = [1 2];
fssSubmesh = [1 2];
ipol = 1;

for igroup=1:fssSolver(fssIndexes(1)).xsLib.ngroups
	j = 0;
	figure(igroup)
	hold on
	tmp = 0;
	for i=fssIndexes
		for k=1:fssSubmesh(i)
			j = j + 1;
			tmpdata(1:length(fssSolver(i).solution.angflux(1,1,ipol,:,k))) = ...
				fssSolver(i).solution.angflux(1,igroup,ipol,:,k);
			plot(celledges,tmpdata,'linestyle',linS{j},'linewidth',2);
			tmp = max(tmp,max(tmpdata));
		end
	end
	ax = gca;
	ax.XAxis.TickValues = xgrid;
	ax.XAxis.MinorTickValues = fssSolver(1).mesh.xsedges;
	ax.XTickLabelRotation = 45;
	ax.FontSize = 11;
	xlabel('Position (cm)');
	ylabel('Angular Flux');
	grid on;
	grid minor;
end