close all;
ncases = length(solutions(1,:));

for i=1:ncases+2
	keffs(2,i) = solutions(2,1).keff(1);
	keffs(3,i) = solutions(3,1).keff(1);
	if i == 1
		vfrac(i) = 0.0;
		keffs(1,i) = solutions(2,1).keff(1);
		keffs(4:7,i) = keffs(1,i);
	elseif i == ncases+2
		vfrac(i) = 100.0;
		keffs(1,i) = solutions(3,i-2).keff(1);
		keffs(4:7,i) = keffs(1,i);
	else
		vfrac(i) = (i-1)/(ncases+1)*100;
		keffs(1,i) = solutions(1,i-1).keff(1);
		keffs(4,i) = solutions(4,i-1).keff(1);
		keffs(5,i) = solutions(5,i-1).keff(1);
		keffs(6,i) = solutions(6,i-1).keff(1);
		keffs(7,i) = solutions(7,i-1).keff(1);
	end
end



isol = 5;
xgrid = 0.5*(fssSolver(1).mesh.fsredges(2:end)+fssSolver(1).mesh.fsredges(1:end-1));
for igroup=1:7
    figure(igroup)
    plot(xgrid,solutions(1,isol).scalflux(igroup,:,1),'k',xgrid,solutions(4,isol).scalflux(igroup,:,1),'r-s',...
        xgrid,solutions(5,isol).scalflux(igroup,:,1),'b:*',xgrid,solutions(6,isol).scalflux(igroup,:,1),'g-.+',...
        xgrid,solutions(7,isol).scalflux(igroup,:,1),'k--o','linewidth',2)
    title(sprintf('Scalar Flux at 50%% Rod Withdrawal, Group %i',igroup))
    legend('Volume Homog.','Subray w/o Recomb.','Subray w/ Recomb. - 0','Subray w/ Recomb. - 1',...
        'Subray w/ Recomb. - 2')
    ax = gca;
    ax.XAxis.TickValues = min(fssSolver(1).mesh.xsedges):input.pitch:max(fssSolver(1).mesh.xsedges);
    ax.XAxis.MinorTickValues = fssSolver(1).mesh.xsedges;
    ax.XTickLabelRotation = 45;
    ax.FontSize = 28;
    xlabel('Position (cm)');
    xlim([0,21.42]);
    ylabel('Scalar Flux');
    grid on;
end
figure(8)
plot(vfrac,keffs(1,:),'k',vfrac,keffs(4,:),'r-s',vfrac,keffs(5,:),'b:*',...
    vfrac,keffs(6,:),'g-.+',vfrac,keffs(7,:),'k--o','linewidth',2)
    legend('Volume Homog.','Subray w/o Recomb.','Subray w/ Recomb. - 0','Subray w/ Recomb. - 1',...
        'Subray w/ Recomb. - 2')
title('K-eff vs. Control Rod Withdrawal');
ax = gca;
ax.FontSize = 28;
xlabel('Control Rod Percent Withdrawn');
ylabel('k-eff');