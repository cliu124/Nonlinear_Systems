function out = analyt_prod(eigL_,eigR)
if size(eigR,2) == 4 || size(eigR,2) == 9 %rotation
	if size(eigL_,2)==3
	inds = [1 2;
	        2 3];
	indA = [1 2;
	        3 4];
	else %3D
	inds = [1 2 4;
	        2 3 5;
	        4 5 6];
	indA = [1 2 3;
	        4 5 6;
	        7 8 9];
	end;
	indB = indA';
	out = zeros(size(eigL_));
	for ii=1:size(inds,1)
	for jj=1:size(inds,1)
	for mm=1:size(inds,1)
	for nn=1:size(inds,1)
		if ii<nn
			continue
		end;
		out(:,inds(ii,nn)) = out(:,inds(ii,nn)) + eigR(:,indB(ii,jj)).*eigL_(:,inds(jj,mm)).*eigR(:,indA(mm,nn));
	end;
	end;
	end;
	end;
else %size(eigR,2) == 3 || size(eigR,2) == 6, coord transform
	if size(eigR,2)==3
	inds = [1 2;
	        2 3];
	else %3D
	inds = [1 2 4;
	        2 3 5;
	        4 5 6];
	end;
	out = zeros(size(eigL_));
	for ii=1:size(inds,1)
	for jj=1:size(inds,1)
		out(:,ii) = out(:,ii) + eigL_(:,jj).*eigR(:,inds(ii,jj));
	end;
	end;
end;
