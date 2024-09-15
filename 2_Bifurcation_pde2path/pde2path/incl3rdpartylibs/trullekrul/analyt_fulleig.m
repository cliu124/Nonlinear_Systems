function out = analyt_fulleig(eigL)
zeron = zeros(size(eigL,1),1);
if size(eigL,2) == 2
	out = [eigL(:,1) zeron eigL(:,2)];
else %3D
	out = [eigL(:,1) zeron eigL(:,2) zeron zeron eigL(:,3)];
end;   