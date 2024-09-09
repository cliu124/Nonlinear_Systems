function [edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri)
edga = [tri(:,2) tri(:,1); tri(:,3) tri(:,2); tri(:,1) tri(:,3)];
edga2tri = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]';

[edg,I] = sort(edga,2); %biggest in last column
[edg,Is] = sortrows(edg);
tris = edga2tri(Is);
I = I(Is,1)==1; %true for non flipped
d =[true; or(edg(1:end-1,1)~=edg(2:end,1), ...
            edg(1:end-1,2)~=edg(2:end,2))]; Nd = sum(d);
ad = not(d); Nad = sum(ad); 
d2 = [ad(2:end); false]; d2 = d2(d);
edg = edg(d,:);

% gen edg2tri
tris2 = zeros(Nd,1); tris2(d2) = tris(ad);
edg2tri = [tris(d) tris2];

%flip edg2tri according to edg permutation
Iflp=and(d2,I(d)); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);

%gen edga2edg
edga2edg = zeros(size(edga,1),1);
edga2edg(Is(d))  = 1:size(edg,1);
edga2edg(Is(ad)) = find(d2);

% gen tri2edga tri2edg
tri2edga = reshape([1:size(tri,1)*3]',size(tri,1),3);
tri2edg = edga2edg(tri2edga);

if nargout <= 3
	return;
end;

% gen nd2edg
nd2edg = inv_table(edg);

if nargout == 8
	return;
end;
% gen nd2tri
nd2tri = inv_table(tri);



% %alt gen nd2tri
% triR = reshape(tri',3*size(tri,1),1); nds=max(tri(:));
% ndC = reshape(repmat((1:size(tri,1)),3,1),3*size(tri,1),1);
% A = sparse(triR,ndC,ndC);
% nN = full(sum(A~=0,2))';
% maxtri = max(nN);
% [nds,R,tris] = find(A);
% C = repmat((1:max(nN))',1,numel(nN)); CnN = repmat(nN,maxtri,1); I = C <= CnN; 
% nd2tri = zeros(nds,maxtri); nd2tri(nds+(C(I)-1)*numel(nN)) = tris;

%alt gen edg
%nds = max(tri(:));
%A = sparse(edga(:,1),edga(:,2),edga2tri);
%[R,C,I] = find(A);
%[R_,C_,I_] = find((A'~=0) + (A~=0));
%edg = [R(R>C) C(R>C)];

