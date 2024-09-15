function out = rpval2M(R,vals)
%C    = cell2mat(accumarray(R,1   ,[],@(x){cumsum(x)})); nMax = max(C);
%vals = cell2mat(accumarray(R,vals,[],@(x){x}));
%R    = cell2mat(accumarray(R,R   ,[],@(x){x}));
%out = zeros(Rmax,nMax); out(R+(C-1)*Rmax) = vals;

%faster alternative using sort()
if not(issorted(R,'rows'))
	[R,I] = sort(R); %[size(vals) size(I)]
	vals = vals(I);
end;
nN = diff([0 find([R(1:end-1)'~=R(2:end)' true])]); nMax = max(nN); Rmax = max(R);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN,nMax,1); I = C <= CnN; 
out = zeros(Rmax,nMax); out(R+(reshape(C(I),size(R))-1)*Rmax) = vals;