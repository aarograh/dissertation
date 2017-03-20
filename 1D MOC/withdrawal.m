TenPin7group_mix_rodded_unrodded
percent = 1;
ncase = 0;
for i=percent:percent:100-percent
	ncase = ncase + 1;
	if ncase == 1
		getref = true;
	else
		getref = false;
	end
	input.mixvols = [i/100, (100-i)/100];
	input.subray = false;
	compare_subrayES;
	solutions(1,ncase) = fssSolver(1).solution;
	solutions(2,ncase) = fssSolver(2).solution;
	solutions(3,ncase) = fssSolver(3).solution;
	solutions(4,ncase) = fssSolver(4).solution;
	close all;
end

save workspaces/TenPin_withdrawal.mat;