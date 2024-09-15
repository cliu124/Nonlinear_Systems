function [ngh,ngh2] = bks_nd2tri2ndO(tri,nd2tri,Ibnd,edg)
% SLOWER ALTERNATIVE FOR ORDERING NGH
tris = size(tri,1);
nds = size(nd2tri,1);
bnds = sum(Ibnd);
Inds = not(Ibnd);
maxnghI = max(sum(nd2tri > 0,2));
maxnghb = max(sum(nd2tri > 0,2)+Ibnd); maxngh = max(maxnghb,maxnghI);
tri_ = [tri; repmat(0,1,size(tri,2))];
trin = nd2tri; trin(trin == 0) = tris+1;
trin = [trin repmat(tris+1,nds,maxnghb-maxnghI)]; %pad tris+1 when maxnghb>maxnghI
triN = sum(trin~=tris+1,2)+Ibnd;
trind = [tri_(trin',1) tri_(trin',2) tri_(trin',3)];
empty = (trin==tris+1)'; 
if size(tri,2) == 4
	trind = [trind tri_(trin',4)];
end;
if nargin==4
cnnds1 = reshape(repmat(edg(:,1)',maxngh*4,1),4,nds*maxngh)';
cnnds2 = reshape(repmat(edg(:,2)',maxngh*4,1),4,nds*maxngh)';
I1 = cnnds1 ~= trind; I1(empty(:),1) = 0;
I2 = cnnds2 ~= trind; I2(empty(:),2) = 0;
I = and(I1,I2);
trind = trind'; trind = reshape(trind(I')',2,size(trind,2))';
I = any([and(I1(:,2) == 0,I2(:,1) == 0) and(I1(:,3) == 0,I2(:,2) == 0) and(I1(:,4) == 0,I2(:,3) == 0) and(I1(:,1) == 0,I2(:,3) == 0) and(I1(:,2) == 0,I2(:,4) == 0) and(I1(:,4) == 0,I2(:,1) == 0)],2);
else %not swap
dimp = size(tri,2);
cnnds = reshape(repmat(1:nds,maxngh*dimp,1),dimp,nds*maxngh)';
I = cnnds ~= trind; I(empty(:),1) = 0;
trind = trind'; trind = reshape(trind(I')',dimp-1,size(trind,2))';
if size(tri,2) == 3
I = I(:,2) == 1;
else %3D
I = or(I(:,2) == 0, I(:,4) == 0);
end;
end;

trind_ = trind(I,1); trind(I,1) = trind(I,2); trind(I,2) = trind_;
if nargout == 2 %size(tri,2) == 4 && nargin~=4 %
ngh=trind; 
ngh2 = reshape(ngh',3,maxngh,nds);
%ngh = reshape(ngh',3,maxngh,nds); ngh2 = mk_tcrcls(ngh(:,:,Ibnd)); return;
[trind,nds,maxngh,triN] = mk_tcrcls(ngh2(:,:,Ibnd));
ngh = zeros(nds,maxngh);
ngh(:,[1 2]) = trind(1:maxngh:end,:);
else
ngh = zeros(nds,maxngh);
ngh(Inds,[1 2]) = trind(reshape([Inds'; false(maxngh-1,nds)],nds*maxngh,1),:);
%find starting boundary node
Ibnd_ = repmat(Ibnd',maxngh,1);
trind_ = reshape(trind(Ibnd_(:),:)',2*maxngh,bnds);
truss_ = repmat([true;false],maxngh,bnds); C = repmat(1:bnds,2*maxngh,1);
[trind_,R] = sort(trind_); truss=zeros(size(R)); truss(:) = truss_(R(:)+(C(:)-1)*2*maxngh);
I = and([trind_(1,:)~=trind_(2,:); and(trind_(2:end-1,:)~=trind_(3:end,:),trind_(2:end-1,:)~=trind_(1:end-2,:)); trind_(end,:)~=trind_(end-1,:)],truss);
thmsk = zeros(2*maxngh,1); thmsk(1:2:2*maxngh) = 1:maxngh;
I_ = false(maxngh,bnds); I_(thmsk(R(I))+(0:bnds-1)'*maxngh) = true;
trind_ = trind(Ibnd_(:),:);
ngh(Ibnd,[1 2]) = trind_(I_(:),:);
end;
R=(1:nds)';

for i=3:maxngh
    C = min(i-1,triN-1);
    ndsrch = repmat(ngh(R+(C-1)*nds)',maxngh,1);
    ngh(R,i) = trind(ndsrch(:) == trind(:,1),2).*(triN>=i);
end;

function [trind,nds,maxngh,Nn] = mk_tcrcls(circles)
%save for_debug3Da.mat;
nds = size(circles,3);
nds_ = size(circles,2);
nd1 = reshape(repmat(1:size(circles,3),3*size(circles,2),1),size(circles));
nd2 = reshape(repmat((1:3)',1,size(circles,2)*size(circles,3)),size(circles));
circles_ = reshape(circles,3*size(circles,2),size(circles,3));
edga2 = circles_(reshape([(2:3:3*nds_); (3:3:3*nds_); (1:3:3*nds_)],3*nds_,1),:);
[edga,I2] = sort([circles_(:) edga2(:)],2);
[edga,I] = sortrows([reshape(nd1,3*nds*nds_,1) edga]);
%d = and(edga(:,2)~=0,any([edga(1,:) ~= edga(2,:); and(edga(2:end-1,:) ~= edga(3:end,:),edga(2:end-1,:) ~= edga(1:end-2,:)); edga(end,:) ~= edga(end-1,:)],2));
d = and(edga(:,2)~=0,[any(edga(1,:) ~= edga(2,:),2); and(any(edga(2:end-1,:) ~= edga(3:end,:),2),any(edga(2:end-1,:) ~= edga(1:end-2,:),2)); any(edga(end,:) ~= edga(end-1,:),2)]);
%save for_debug3D.mat; error('just stop3D');
edga = edga(d,:);
Iflp = I2(I(d),1) == 2;
edga(Iflp,:) = edga(Iflp,[1 3 2]);
Nn = diff([0; find([edga(1:end-1,1)~=edga(2:end,1); true])]); maxngh = max(Nn);
R = ones(size(edga,1),1); R(1+cumsum(Nn(1:end-1))) = 1+maxngh-Nn(1:end-1);
R = cumsum(R);
trind = zeros(nds*maxngh,2);
trind(R,:) = edga(:,[2 3]);

%function ngh = ngh_ndsO(edg,tri,nd2edg,nd2tri,edg2tri)
%% FASTEST ALTERNATIVE FOR ORDERED NGH
%tris = size(tri,1);
%nds = size(nd2tri,1);
%edgs = size(edg,1);
%bndedg = sum(edg2tri>0,2)==1;
%bndnd = edg(bndedg,:); bndnd = sort(bndnd(:)); bndnd = bndnd(find([diff(bndnd)>0; 1]));
%bndrnt = zeros(size(nd2tri,1),1); bndrnt(bndnd) = 1;
%maxnghI = max(sum(nd2tri > 0,2));
%maxnghb = max(sum(nd2tri > 0,2)+bndrnt); maxngh = max(maxnghb,maxnghI);
%Ibnd = find(bndrnt);
%Inds = find(not(bndrnt));
%Iall = 1:nds;
%ncorner = sum(nd2tri(Ibnd,:)>0,2)>1;
%tri_ = [tri; 0 0 0];
%trin = [nd2tri(Iall,:) zeros(nds,max(0,maxnghb-maxnghI))]; trin(trin == 0) = tris+1;
%triN = sum(trin~=tris+1,2);
%trind = [tri_(trin',1) tri_(trin',2) tri_(trin',3)];
%cnnds = reshape(repmat(Iall,maxngh*3,1),3,nds*maxngh)';
%I = cnnds ~= trind;
%empty = (trin==tris+1)'; 
%I(empty(:),1) = 0;
%trind = trind'; trind = reshape(trind(I')',2,size(trind,2))';
%I = I(:,2) == 1;
%% swapping trind(I,1) and trind(I,2)
%trind_ = trind(I,1); trind(I,1) = trind(I,2); trind(I,2) = trind_;
%ngh = zeros(nds,maxngh);
%ngh(Inds,1) = trind(1+(Inds-1)*maxngh,1);
%nd2edg_ = nd2edg(Ibnd,:)'; nd2edg_(nd2edg_==0) = edgs + 1; bndedg_ = [bndedg; false]; nd2edg_ = nd2edg_(bndedg_(nd2edg_)); %cmp = ismember(nd2edg_,find(bndedg)); nd2edg_ = nd2edg_(cmp);
%bndndsp1 = edg(nd2edg_(1:2:end),:)'; cmp = repmat(Ibnd',2,1)~=bndndsp1; bndndsp1 = bndndsp1(cmp);
%bndndsp2 = edg(nd2edg_(2:2:end),:)'; cmp = repmat(Ibnd',2,1)~=bndndsp2; bndndsp2 = bndndsp2(cmp);
%ndsrch = repmat(bndndsp2',maxngh,1);
%I = repmat((1:maxngh)',1,numel(Ibnd))+repmat((Ibnd-1)*maxngh,1,maxngh)';
%I2 = find(trind(I(:),1) == ndsrch(:));
%ngh(Ibnd,1) = bndndsp1;
%ngh(Ibnd((I2-mod(I2,maxngh))/maxngh+1),1) = bndndsp2((I2-mod(I2,maxngh))/maxngh+1);
%triN(Ibnd) = triN(Ibnd) + 1; triN(Ibnd(not(ncorner))) = 2;
%triNsz = size(trind,1)/nds;
%R = Iall';
%for i=2:maxngh
    %C = min(i-1,triN-1);
    %ndsrch = repmat(ngh(R+(C-1)*nds)',triNsz,1); %maxngh,1);
    %ngh(R,i) = trind(ndsrch(:) == trind(:,1),2).*(triN>=i);
%end;

%[R,Cc] = find(ngh~=0); nN = sum(ngh~=0,2); msk1 = [0:size(ngh,2)-1]'; msk2 = [2:size(ngh,2)+1]';
%Cl = msk1(Cc); Cl(Cl==0) = nN(R(Cl==0));
%Cr = msk2(Cc); Cr(Cr==nN(R)+1) = 1;
%newtri = [[edg(R,1) ngh(R+(Cl-1)*size(ngh,1)) ngh(R+(Cc-1)*size(ngh,1)) ngh(R+(Cr-1)*size(ngh,1))]; [edg(R,2) ngh(R+(Cc-1)*size(ngh,1)) ngh(R+(Cl-1)*size(ngh,1)) ngh(R+(Cr-1)*size(ngh,1))]];
%[tmp,area] = elem_inv(newtri,xy);
%result = false([size(ngh,1) size(ngh,2)*2]);
%result(R+(2*Cc-2)*size(ngh,1)) = tmp(1:numel(R));
%result(R+(2*Cc-1)*size(ngh,1)) = tmp(numel(R)+1:end); sum(not(any(result,2)))

