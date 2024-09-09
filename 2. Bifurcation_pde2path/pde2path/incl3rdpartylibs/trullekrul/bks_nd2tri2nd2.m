function ngh = bks_nd2tri2nd2(tri,nd2tri,Ind)
%SIMPLE ALTERNATIVE WITHOUT ORDERING
if nargin == 2
	Ind = 1:size(nd2tri,1);
end;
nd2tri = nd2tri';
I = nd2tri~=0; nd2tri_ = nd2tri(I);
ngh1 = zeros(size(nd2tri)); ngh1(I) = tri(nd2tri_,1);
ngh2 = zeros(size(nd2tri)); ngh2(I) = tri(nd2tri_,2);
ngh3 = zeros(size(nd2tri)); ngh3(I) = tri(nd2tri_,3);
ngh = [ngh1; ngh2; ngh3];
if size(tri,2) == 4
 ngh4 = zeros(size(nd2tri)); ngh4(I) = tri(nd2tri_,4);
 ngh = [ngh; ngh4];
end;
ngh(ngh==repmat(Ind,size(ngh,1),1)) = 0;
ngh = sort(ngh);
ngh([false(1,size(ngh,2)); ngh(1:end-1,:)==ngh(2:end,:)]) = 0;
nN = sum(ngh~=0);
nds = ngh(ngh~=0);
ngh = zeros(max(nN),size(ngh,2));
ngh(repmat((1:max(nN))',1,size(ngh,2)) <= repmat(nN,max(nN),1)) = nds;
ngh = ngh';