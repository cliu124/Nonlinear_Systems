function [xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_rm_edg(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
% GENERATE BOOKS
ndone = 0;
if size(xy,2) == 2
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg,nd2tri] = bks_all(tri);
[nghOa,ngh,deltri] = bks_edg2nd2tri2nd2edg(edg,tri,nd2edg,nd2tri,edga,edga2edg,tri2edga,options);
	else
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg,nd2tri] = bks_all3D(tri);
[tmp,edg1,edg2] = geom_crnds(bndmesh);
edga = faca;
[nghOa,ngh,deltri] = bks_edg2nd2tri2nd2edg(edg,tri,nd2edg,nd2tri,faca,faca2fac,tri2faca,options,nd2fac);
end;

edg2clr = bks_clr(ngh,options); 

clrs = max(edg2clr);
bndedg_ = repmat(size(xy,2),size(edg,1),1); %zeros(size(edg,1),1); 
bndnds_ = repmat(size(xy,2),size(xy,1),1); %zeros(size(xy,1),1); 
if size(xy,2) == 2
bndedg2edg = bks_bndedg2edg(edg,nd2edg,bndmesh.edg);
bndedg_(bndedg2edg) = 1;
bndnds_(bndmesh.edg(:)) = 1;
bndnds_(bndmesh.crnds) = 0;
else %3D
bndedg_(bks_bndedg2edg(edg,nd2edg,edg2)) = 2;
bndnds_(bndmesh.fac(:)) = 2;
bndedg_(bks_bndedg2edg(edg,nd2edg,edg1)) = 1;
bndnds_(edg1(:)) = 1;
bndnds_(bndmesh.crnds) = 0;
end;

%bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = 1;
if options.fastRM
	[rbadedg,ctrs] = find_badedg(edg,xy,Nmetric,options); %FIND BAD EDGES
else
	rbadedg = 1:size(edg,1); [tmp,ctrs] = elem_qual(edg,xy,Nmetric,options); %we ought to implement an option to avoid tmp computation
end;
%if options.qualM == 9
%	Nmetric_ = Nmetric;
%	Nmetric = metric_iso(xy, 1);
%end;
%do not contract edge lying in D object if both nodes are part of a D-1 object
Igood = max([bndnds_(edg(rbadedg,1)) bndnds_(edg(rbadedg,2))],[],2)>=bndedg_(rbadedg);
rbadedg = rbadedg(Igood); ctrs = ctrs(Igood,:);
%save for_debug.mat;
ctrs_ = zeros(size(edg,1),size(xy,2)); ctrs_(rbadedg,:) = ctrs;

%update colours to reflect bad edges
rbadedg_ = false(size(edg,1),1); rbadedg_(rbadedg) = true; edg2clr(not(rbadedg_)) = 0;

delnds = [];
deledgs = [];
ebadnds_ = [];
deltri_ = [];
newtri = [];
newQ = [];
newxy_ = [];
newNmetric_ = [];
collapsends = 1:size(xy,1);

for i=1:clrs
    badedg = find(edg2clr==i);
    if numel(badedg) == 0 %nothing to do here
        continue;
    end;
    [newtri_,newxy,badedg,newNmetric,Nqual] = mesh_bnd_interior(nghOa(badedg,:),edga,edg(badedg,:),xy,Nmetric,badedg,triQ,deltri(badedg,:),bndnds_,ctrs_(badedg,:),geomfunc,options);  
    if numel(badedg) == 0 %we could not improve the situation
        continue;
    else
        ndone = ndone + numel(badedg);
    end;
    newnds = repmat((1:size(newxy,1))',1,2)+size(xy,1)+size(newxy_,1);
    %newnds = repmat(1:size(newxy,1),2,1)+size(xy,1)+size(newxy_,1);
    badnds = edg(badedg,:); badnds = badnds(:);
    collapsends(badnds) = newnds(:);
    ebadnds = extra_badnds(deltri(badedg,:),nghOa(badedg,:),edga,tri,edg(badedg,:),options);
    delnds = [delnds; badnds];
    deledgs = [deledgs; badedg];
    ebadnds_ = [ebadnds_; ebadnds];
    newtri_(:,2) = newtri_(:,2) + size(newxy_,1);
    newxy_ = [newxy_; newxy];
    newNmetric_ = [newNmetric_; newNmetric];
    newtri = [newtri; newtri_];
    newQ = [newQ; Nqual];
    if i < clrs %update without running gen_books
        clrlssedg = ngh(badedg,:); clrlssedg = clrlssedg(clrlssedg~=0);
        edg2clr(clrlssedg(:)) = 0; 
    end;
end;
%if options.qualM == 9
%	Nmetric = Nmetric_;
%end;
if ndone ~= 0
ndone = ndone/size(xy,1);
I = true(size(xy,1),1); I([delnds; ebadnds_]) = false; xy = [xy(I,:); newxy_]; Nmetric = [Nmetric(I,:); newNmetric_];
old2new = mvndtri([I; true(size(newxy_,1),1)]); %numel(I)+size(newxy_,1),
deltri_ = deltri(deledgs,:); deltri_ = deltri_(deltri_~=0);
if options.mntn_bks
bks.mvnd(tri(deltri_,:)) = true;
bks.rmnd(tri(deltri_,:)) = true;
bks.mvnd = [bks.mvnd(I); true(size(newxy_,1),1)];
bks.rmnd = [bks.rmnd(I); true(size(newxy_,1),1)];
end;

I = true(size(tri,1),1); I(deltri_) = false;
tri = [tri(I,:); newtri];
triQ = [triQ(I,:); newQ];
tri = old2new(tri);
%update boundary mesh
%I_ = zeros(size(I)); I_(delnds) = %reshape(repmat(1:size(newxy_,1),2,1),size(delnds));
%I_ = I_(bndmesh.crnds); bndmesh.crnds(I_~=0) = I_(I_~=0)+size(I,1);
bndmesh.crnds = collapsends(bndmesh.crnds);
bndmesh.crnds = old2new(bndmesh.crnds);
if size(xy,2) == 3
	bndmesh.fac = collapsends(bndmesh.fac);
	bndmesh = fix_bndmesh(bndmesh,ebadnds_);
	%bndmesh = rmfield(bndmesh,'edg');
	bndmesh.fac = old2new(bndmesh.fac);
	bndmesh.fac = sort(bndmesh.fac,2);
else
	I = true(size(edg,1),1); I(deledgs) = false;
	keepbndedg = I(bndedg2edg);
	bndmesh.edg = collapsends(bndmesh.edg(keepbndedg,:)); bndmesh.IDs = bndmesh.IDs(keepbndedg); 
	bndmesh.edg = old2new(bndmesh.edg);
	bndmesh.edg = sort(bndmesh.edg,2);
end;
if options.debug
    sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;


function [newtri,newxy,badedg,newNmetric,newqual] = mesh_bnd_interior(nghOa_,edga,edgbad,xy,Nmetric,badedg,triQ,deltri,bndnds_,ctrs,geomfunc,options)
newxy = calc_newxy(xy,edgbad,badedg,bndnds_,ctrs,geomfunc);
newNmetric = metric_avg(edgbad,Nmetric,options); 
xy_ = [xy; newxy];
Nmetric_ = [Nmetric; newNmetric];
newxys = repmat((1:numel(badedg))'+size(xy,1),1,size(nghOa_,2));
I = nghOa_~=0; Iedga = nghOa_(I);
newtri = [edga(Iedga,1) reshape(newxys(I),sum(I(:)),1) edga(Iedga,2)];
if size(edga,2) == 3
	newtri = [newtri(:,[3 2 1]) edga(Iedga,3)];
end;
%[I_,z] = elem_inv(newtri,xy_); Ibad = z<options.minA; %mean(sign(z))
Nqual = elem_qual(newtri,xy_,Nmetric_,options,1);
R = repmat((1:size(nghOa_,1))',1,size(nghOa_,2)); C = repmat(1:size(nghOa_,2),size(nghOa_,1),1); R = R(I); C=C(I);
%Iinv = false(size(nghOa_)); Iinv(R+(C-1)*size(Iinv,1)) = Ibad;
Nquality= inf(size(nghOa_)); Nquality(R+(C-1)*size(nghOa_,1)) = Nqual;
%notinv = not(any(Iinv,2));
worstN = min(Nquality,[],2);
if options.consRM
	deltri_ = deltri; deltri_(deltri_==0) = size(triQ,1)+1; triQ_ = [triQ; Inf];
	triQ_ = reshape(triQ_(deltri_),size(nghOa_,1),size(deltri_,2));
	worstO = min(triQ_,[],2);
	goodnotinv = worstN > worstO;
	%notinv(notinv) = goodnotinv; 
else
	goodnotinv = worstN > 0;
end;	
badedg = badedg(goodnotinv);
newNmetric = newNmetric(goodnotinv,:);
newqual = Nquality(goodnotinv,:);
newxy = newxy(goodnotinv,:);
%restate newtri (speedier alternative might exist)
nghOa_=nghOa_(goodnotinv,:);
In = nghOa_~=0; Iedga = nghOa_(In);
newxys = repmat((1:numel(badedg))'+size(xy,1),1,size(nghOa_,2));
newtri = [edga(Iedga,1) reshape(newxys(In),sum(In(:)),1) edga(Iedga,2)];
if size(edga,2) == 3
	newtri = [newtri(:,[3 2 1]) edga(Iedga,3)];
end;
newqual = reshape(newqual(not(isinf(newqual))),sum(In(:)),1);


function newxy = calc_newxy(xy,edgbad,badedg,bndnds_,ctrs,geomfunc)
%identify edges related to boundaries
I = repmat(2,numel(badedg),1); %interior edge is the default
Idim = [bndnds_(edgbad(:,1)) bndnds_(edgbad(:,2))];
I(and(Idim(:,1)==Idim(:,2),Idim(:,1)~=size(xy,2))) = 1; %collapse to center
I(Idim(:,1)~=Idim(:,2)) = 0; %collapse to node part of lower dimensional
newxy = zeros(numel(badedg),size(xy,2));
Isimp = or(I==2,I==1);
newxy(Isimp ,:) = ctrs(Isimp,:);
if numel(geomfunc)~=0 && any(I==1)
	newxy_ = newxy(I==1,:);
	dist = geomfunc(newxy_ );
	newxy(I==1,:) = newxy_ - dist(:,[2 3]).*repmat(dist(:,1),1,2);
end;
if any(I==0)
Icomp = I==0;
edg_ = edgbad(Icomp,:)';
newxy(Icomp,:) = xy(edg_(Idim(Icomp,:)' < repmat(max(Idim(Icomp,:)'),2,1)),:);
end;

function badnds = extra_badnds(deltri,nghOa_,edga,tri,badedga,options)
nds = max(edga(:)); Niedg = size(deltri,1);
tri = tri';
iedg = (1:size(deltri,1))'; R = repmat(iedg,1,size(deltri,2));
R = R(deltri~=0); deltri = deltri(deltri~=0); R = reshape(R,numel(R),1); deltri = reshape(deltri,numel(deltri),1);
nds1 = tri(:,deltri); R3 = repmat(R',size(tri,1),1); 
if options.ismatlab %numel(strfind(version,'R2010a'))~=0 %|| options.nosparse
%take out doubles
[nds1R3,I] = sortrows([nds1(:) R3(:)]);
I = [true; or(nds1R3(2:end,1)~=nds1R3(1:end-1,1),nds1R3(2:end,2)~=nds1R3(1:end-1,2))];
nds1R3 = nds1R3(I,:);
end;
if options.nosparse
	tbl1 = [nds1(:) R3(:)]; %tbl1 = [nds1R3(:,1) nds1R3(:,2)];
else
if options.ismatlab  %numel(strfind(version,'R2010a'))~=0
A = sparse(nds1R3(:,1),nds1R3(:,2),true,nds,numel(iedg)); 
else
A = sparse(nds1(:),R3(:),true,nds,numel(iedg)); 
end;
end;
R = repmat(iedg,1,size(nghOa_,2));
R = R(nghOa_~=0); nghOa_ = nghOa_(nghOa_~=0); R = reshape(R,numel(R),1); nghOa_ = reshape(nghOa_,numel(nghOa_),1);
nds2 = edga(nghOa_,:);  R2 = repmat(R,size(edga,2),1);
if options.ismatlab %numel(strfind(version,'R2010a'))~=0 %|| options.nosparse
%take out doubles
nds2 = [nds2(:); badedga(:,1); badedga(:,2)]; R2=[R2;iedg;iedg];
[nds2R2,I] = sortrows([nds2 R2]);
I = [true; or(nds2R2(2:end,1)~=nds2R2(1:end-1,1),nds2R2(2:end,2)~=nds2R2(1:end-1,2))];
nds2R2 = nds2R2(I,:);
end;
if options.nosparse
	tbl2 = [[nds2(:); badedga(:,1); badedga(:,2)] [R2;iedg;iedg]];%tbl2 = [nds2R2(:,1) nds2R2(:,2)]; 
	nvtbl2 = inv_table(tbl2(:,2));
	nvtbl_ = nvtbl2(tbl1(:,2),:); I = nvtbl_~=0; nvtbl_(I) = tbl2(nvtbl_(I),1);
	cmp = nvtbl_ == repmat(tbl1(:,1),1,size(nvtbl2,2));
	badnds = tbl1(not(any(cmp,2)),1);
else
if options.ismatlab  %numel(strfind(version,'R2010a'))~=0
A2 = sparse(nds2R2(:,1),nds2R2(:,2),true,nds,Niedg);
else
A2 = sparse([nds2(:); badedga(:,1); badedga(:,2)],[R2;iedg;iedg],true,nds,Niedg);
end;
%A2 = sparse([nds2(:); badedga(:,1); badedga(:,2)],[R;R;iedg;iedg],true,nds,Niedg); [C2,R2] = find(A2);
%A2=sparse(C2,R2,1,nds,Niedg);
A = A-A2.*A; [badnds,R,I] = find(A); %badnds = badnds(I~=0);
end; 
if numel(badnds) > 1
	badnds = sort(badnds); badnds = badnds([true; 			         badnds(2:end)~=badnds(1:end-1)]);
	badnds = reshape(badnds,numel(badnds),1);
elseif numel(badnds) == 0
	badnds = []; %avoids 0x1 size
end;

function tri = mvndtri(goodnds) %reduces node numbers
tri = zeros(size(goodnds)); tri(goodnds) = 1:nnz(goodnds);

%function tri = mvndtri(maxtri,goodnds) %reduces node numbers
%badnds = badnds(badnds<=maxtri);
%if maxtri >= min(badnds) 
%    tri = 1:maxtri;
%    move = zeros(1,max(tri)); move(badnds) = - 1; move = cumsum(move);
%    tri = 1:max(tri(:)); tri = tri + move;
%end;

function [badedg,ctrs] = find_badedg(edg,xy,Nmetric,options)
[edgnowl,ctrs] = elem_qual(edg,xy,Nmetric,options);
I = options.Llow > edgnowl;
badedg = find(I);
ctrs = ctrs(I,:);

function bndmesh = fix_bndmesh(bndmesh,ebadnds)
I = any([bndmesh.fac(:,1) == bndmesh.fac(:,2) bndmesh.fac(:,1) == bndmesh.fac(:,3) bndmesh.fac(:,2) == bndmesh.fac(:,3)],2);
if numel(ebadnds) == 0
	nds = max(bndmesh.fac(:));
else
	nds = max(max(ebadnds),max(bndmesh.fac(:)));
end;
ebadnds_ = false(nds,1); ebadnds_(ebadnds) = true;
I2 = any(ebadnds_(bndmesh.fac),2);
I = not(or(I,I2));
bndmesh.fac = bndmesh.fac(I,:);
bndmesh.IDs = bndmesh.IDs(I,:);