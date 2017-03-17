TenPin7group_mix_rodded_unrodded
percent = 5;
ncase = 0;
for i=percent:percent:100-percent
	ncase = ncase + 1;
	input.mixvols = [i/100, (100-i)/100];
	compare_subrayES;
	solutions(1,ncase) = fssSolver(1).solution;
	solutions(2,ncase) = fssSolver(2).solution;
	solutions(3,ncase) = fssSolver(3).solution;
	solutions(4,ncase) = fssSolver(4).solution;
	clear eSolver;
	close all;
end

save workspaces/TenPin_withdrawal.mat;