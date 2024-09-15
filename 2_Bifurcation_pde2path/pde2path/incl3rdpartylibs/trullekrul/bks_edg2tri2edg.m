function edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg,Iedg)
tris = size(tri2edg,1);
edgs = size(edg2tri,1);
if nargin == 2
	Iedg = (1:edgs)';
end;

switch size(tri2edg,2)
case 3
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
N1 = zeros(edgs,3); N1(I(:,1),:) = tri2edg(edg2tri(I(:,1),1),:);
N2 = zeros(edgs,3); N2(I(:,2),:) = tri2edg(edg2tri(I(:,2),2),:);
edg_ngh = [N1 N2];
edg_ngh(edg_ngh == repmat(Iedg,1,6)) = 0;
edg_ngh = sort(edg_ngh,2); edg_ngh = edg_ngh(:,3:6);
case 4 %3D face, fac2tri2fac
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
N1 = zeros(edgs,4); N1(I(:,1),:) = tri2edg(edg2tri(I(:,1),1),:);
N2 = zeros(edgs,4); N2(I(:,2),:) = tri2edg(edg2tri(I(:,2),2),:);
fac_ngh = [N1 N2];
%fac_ngh = [tri2edg_(edg2tri_(:,1),:) tri2edg_(edg2tri_(:,2),:)];
fac_ngh(fac_ngh == repmat(Iedg,1,8)) = 0;
edg_ngh = sort(fac_ngh,2); edg_ngh = edg_ngh(:,3:8);%3D swapping
case 6 %3D edg
edg2tri = edg2tri';
tri2edg = tri2edg';
I = edg2tri~=0;
edg2tri_ = edg2tri(I);
edg_ngh = zeros(6*size(edg2tri,1),edgs); edg_ngh(repmat(I,6,1)) = tri2edg(:,edg2tri_);
edg_ngh(edg_ngh == repmat(Iedg',size(edg_ngh,1),1)) = 0;
edg_ngh = sort(edg_ngh);
I = and(edg_ngh~=0,[edg_ngh(1:end-1,:)~=edg_ngh(2:end,:); true(1,edgs)]);
nN = sum(I);
%vals = edg_ngh(I);
%R = edg_(I);
%edg_ngh = rpval2M(vals,R);
edg_ngh(not(I)) = 0;
edg_ngh = sort(edg_ngh);
edg_ngh = edg_ngh(end-max(nN)+1:end,:)';
end;
%
%tris = size(tri2edg,1);
%edgs = size(edg2tri,1);
%edg2tri_ = edg2tri; edg2tri_(edg2tri_==0) = tris+1;
%tri2edg_ = [tri2edg; zeros(1,size(tri2edg,2))]; 
%switch size(tri2edg,2)
%case 3
%edg_ngh = [tri2edg_(edg2tri_(:,1),:) tri2edg_(edg2tri_(:,2),:)];
%edg_ = repmat((1:edgs)',1,6); edg_ngh(edg_ngh == edg_) = 0;
%edg_ngh = sort(edg_ngh,2); edg_ngh = edg_ngh(:,3:6);
%case 4 %3D face, fac2tri2fac
%fac_ngh = [tri2edg_(edg2tri_(:,1),:) tri2edg_(edg2tri_(:,2),:)];
%fac_ = repmat((1:edgs)',1,8); fac_ngh(fac_ngh ==fac_) = 0;
%edg_ngh = sort(fac_ngh,2); edg_ngh = edg_ngh(:,3:8);%3D swapping
%case 6 %3D edg
%edg2tri_ = edg2tri_';
%edg_ngh = reshape(tri2edg_(edg2tri_(:),:)',size(tri2edg,2)*size(edg2tri,2),edgs);
%edg_ = repmat((1:edgs),size(edg_ngh,1),1); edg_ngh(edg_ngh == edg_) = 0;
%edg_ngh = sort(edg_ngh);
%I = and(edg_ngh~=0,[edg_ngh(1:end-1,:)~=edg_ngh(2:end,:); true(1,edgs)]);
%nN = sum(I);
%%vals = edg_ngh(I);
%%R = edg_(I);
%%edg_ngh = rpval2M(vals,R);
%edg_ngh(not(I)) = 0;
%edg_ngh = sort(edg_ngh);
%edg_ngh = edg_ngh(end-max(nN)+1:end,:)';
%end;