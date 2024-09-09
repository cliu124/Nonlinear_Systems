function [edga_,tri2edga_,edga2edg_,tri2edg_,edg2tri_] = bks_all_(tri,edg,nd2tri,nd2edg)
edga_ = [tri(:,2) tri(:,1); tri(:,3) tri(:,2); tri(:,1) tri(:,3)];
tri2edga_ = reshape([1:size(tri,1)*3]',size(tri,1),3);
if nargout == 2
	return;
end;
%gen edga2edg (and tri2edg)
edga2tri_ = [1:size(tri,1) 1:size(tri,1) 1:size(tri,1)]';
[edgas,I] = sort(edga_,2);
cmp = nd2edg(edgas(:,1),:)';
edg2 = zeros(size(nd2edg,2),size(edga_,1)); edg2(cmp~=0) = edg(cmp(cmp~=0),2);
I = repmat(edgas(:,2)',size(nd2edg,2),1) == edg2;
edga2edg_ = cmp(I);
tri2edg_ = edga2edg_(tri2edga_);
if nargout == 4
	return;
end;
%gen edg2tri
cmp = nd2tri(edg(:,1),:)';
tri1 = zeros(size(nd2tri,2),size(edg,1)); tri1(cmp~=0) = tri(cmp(cmp~=0),1);
tri2 = zeros(size(nd2tri,2),size(edg,1)); tri2(cmp~=0) = tri(cmp(cmp~=0),2);
tri3 = zeros(size(nd2tri,2),size(edg,1)); tri3(cmp~=0) = tri(cmp(cmp~=0),3);
I1a = repmat(edg(:,2)',size(nd2tri,2),1) == tri1;
I2a = repmat(edg(:,2)',size(nd2tri,2),1) == tri2;
I3a = repmat(edg(:,2)',size(nd2tri,2),1) == tri3;
I1b = repmat(edg(:,1)',size(nd2tri,2),1) == tri1;
I2b = repmat(edg(:,1)',size(nd2tri,2),1) == tri2;
I3b = repmat(edg(:,1)',size(nd2tri,2),1) == tri3;
I1 = and(I1b,I2a); I2 = and(I2b,I3a); I3 = and(I3b,I1a);
I4 = and(I1a,I2b); I5 = and(I2a,I3b); I6 = and(I3a,I1b);
edg2tri_ = zeros(size(edg,1),2);
edg2tri_(any(I1),1) = cmp(I1); edg2tri_(any(I4),2) = cmp(I4);
edg2tri_(any(I2),1) = cmp(I2); edg2tri_(any(I5),2) = cmp(I5);
edg2tri_(any(I3),1) = cmp(I3); edg2tri_(any(I6),2) = cmp(I6);
I = edg2tri_(:,1) == 0; edg2tri_(I,:) = edg2tri_(I,[2 1]);

%times = [0 0]; tic; [edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri); times(1) = toc; tic; [edga,tri2edga,edga2edg,tri2edg,edg2tri] = bks_all_(tri,edg,nd2tri,nd2edg); times(2) = toc; times(1)/times(2)