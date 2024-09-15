function [mmI,eigL,eigR] = metric_sqrt(mm,operstr)
[eigL,eigR] = analyt_eig(mm);
if nargin == 1
	eigL = sqrt(eigL);
else
	if strcmp(operstr,'inv')
		eigL = 1./eigL;
	elseif strcmp(operstr,'invsqrt')
		eigL = 1./sqrt(eigL);
	else
		error('invalid 2nd function argument')
	end;
end;
mmI = analyt_prod(analyt_fulleig(eigL),eigR);