function ngh = bks_nd2tri2nd(edg,nd2edg,Ind)
%% SIMPLE 7 LINE ALTERNATIVE WITHOUT ORDERING
if nargin == 2
	Ind = (1:size(nd2edg,1))';
end;
I = nd2edg ~= 0; nd2edg_ = nd2edg(I);
N1 = zeros(size(nd2edg)); N1(I) = edg(nd2edg_,1);
N2 = zeros(size(nd2edg)); N2(I) = edg(nd2edg_,2);
tmp = repmat(Ind,1,size(N1,2));
ngh = zeros(size(nd2edg));
I1 = N1~=tmp; ngh(I1) = N1(I1);
I2 = N2~=tmp; ngh(I2) = N2(I2);