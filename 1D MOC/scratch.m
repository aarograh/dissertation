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

plot(vfrac,keffs(1,:),'k',vfrac,keffs(4,:),'r--s',vfrac,keffs(5,:),'b:*',...
    vfrac,keffs(6,:),'g-.+',vfrac,keffs(7,:),'k-o')
legend('Volume Homogenization','Subray','Subray w/ Recombination - 0','Subray w/ Recombination - 1',...
    'Subray w/ Recombination - 2')