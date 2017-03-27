TenPin7group_mix_rodded_unrodded
% ThreePin7group_mix_rodded_unrodded
percent = 10;
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
    for j=1:length(fssSolver)
        solutions(j,ncase) = fssSolver(j).solution;
    end
	close all;
end

afsdiasdf

% save workspaces/TenPin_withdrawal.mat;