function [tri,triQ,bks,bndmesh,ndone,xy,Nmetric] = adapt_rm_nd(tri,xy,Nmetric,bndmesh,triQ,bks,options)
%if options.debug
	%sanity_check(tri,xy,triQ,Nmetric,options);
%end;
% GENERATE BOOKS
ndone = 0;
nd2tri = inv_table(tri);
if size(xy,2) == 3		 
        [fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2] = bks_all3D(tri);
	nd2bndfac = inv_table(bndmesh.fac);
	nd2bndfac = [nd2bndfac; zeros(size(xy,1)-size(nd2bndfac,1),size(nd2bndfac,2))];
	%nd2edg = inv_table(edg);
	bndnds_ = false(size(xy,1),1); bndnds_(fac(fac2tri(:,2)==0,:)) = true;%bndnds_(bndmesh.fac(:)) = true;
	[ngh3,ngh4] = bks_nd2tri2ndO(tri,nd2tri,bndnds_);
	ngh = bks_nd2tri2nd2(tri,nd2tri); %for colouring
	%bndnds = zeros(size(bndnds_)); bndnds(bndnds_) = find(bndnds_);
	[tmp,edg1,edg2,edg2fac,edg2ID]  = geom_crnds(bndmesh);
%	bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
	bndnds2_ = false(size(xy,1),1); bndnds2_(edg1) = true;
	if isfield(bndmesh,'triID') && size(edg2ID,2) > 2
		bndnds3_ = false(size(xy,1),1); bndnds3_(edg2(edg2ID(:,3)~=0,:)) = true; %belonging to three or more faces
		bndndsI_ = false(size(xy,1),1); bndndsI_(bndmesh.fac(:)) = true;	
		ngh3 = ngh3(not(bndnds3_(bndnds_)),:);
		%bndnds2_(and(bndndsI_,not(bndnds_))) = false; %killing the internal ones
		bndnds_(bndnds3_) = false;
		%bndndsI_ = and(bndndsI_,or(bndnds3_,not(bndnds_)));
		bndndsI_(bndnds_) = false; %belonging to interior face		
	end;
	if size(Nmetric,2) ~= 6
	Tmetric = Nmetric(:,7:end); Nmetric = Nmetric(:,1:6);
	end;
else
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	if isfield(bndmesh,'triID')
		[edg,edg2tri,tri2edg] = bks_all(tri);
		bndnds2_ = bndnds_; bndnds2_(edg(edg2tri(:,2)==0,:)) = false;
		ngh = bks_nd2tri2ndO(tri,nd2tri,and(bndnds_,not(bndnds2_)));
		nd2bndedg = inv_table(bndmesh.edg);
		allIDs = bndmesh.triID(nd2tri(:,1));
	else
		ngh = bks_nd2tri2ndO(tri,nd2tri,bndnds_);
		bndnds2_ = false(size(bndnds_));
	end;
	if size(Nmetric,2) ~= 3
	Tmetric = Nmetric(:,4:end); Nmetric = Nmetric(:,1:3);
	end;
end;
%Do not remove corner nodes
nd2clr = bks_clr(ngh,options); nd2clr(bndmesh.crnds) = 0;
if isfield(options,'TS')
	if size(xy,2) == 3
		nd2clr(bndmesh.fac(bndmesh.IDs==options.TS,:)) = 0;
	else
		nd2clr(bndmesh.edg(bndmesh.IDs==options.TS,:)) = 0;
	end;
end;

if options.fastRM
    if options.fastRM == 2
    if size(xy,2) == 2
    edg = bks_all(tri);
    end;
    badedg = find_badedg(xy,Nmetric,edg,options); %FIND BAD EDGES
    rbadnds = badedg2badnds(badedg,edg,size(xy));
    elseif options.fastRM == 1
    rbadnds = find_badnds(xy,Nmetric,ngh,options);
    else
    rbadnds = (1:size(xy,1))';
    end;
    %update colours to reflect bad nodes
    rbadnds_ = true(size(ngh,1),1); rbadnds_(rbadnds) = false; 
    %rbadnds_(and(not(bndnds_),sum(ngh~=0,2)<=size(tri,2)^2-3)) = false; 
    nd2clr(rbadnds_) = 0;
end;

if options.consRM
  triQtb = inf(size(nd2tri)); I = nd2tri~=0;
  triQtb(I) = triQ(nd2tri(I)); triQtb = min(triQtb,[],2);  
else
  triQtb = [];
end;

newtri = [];
newQ = [];
newIDt = [];
delnds = false(size(nd2tri,1),1); 
ndelnds = false(size(nd2tri,1),1);
newfacs = []; newIDs = []; newtriN2 = []; %for 3D coarsening only
clrs = max(nd2clr);
for i=1:clrs
    badnds_ = nd2clr==i; 
    if options.mntn_bks
    	badnds_ = and(badnds_,bks.rmnd);
    end;
    badnds = find(badnds_);
    if numel(badnds) == 0 %nothing to do here
        continue;
    end;
    if size(xy,2) == 2
    	circles = ngh(badnds,:); 
    	edg_ = [];
    else
    	spheres = ngh4(:,:,badnds_);
    	spheres = spheres(:,:,1:max(sum(spheres(1,:,:)~=0,3)));
    	Ifac = and(bndnds_,badnds_); 
    	if any(Ifac)
    	badnds_fac = find(Ifac);
    	circles = ngh3(badnds_(bndnds_),:);
%    	if numel(geomfunc) ~= 0
%    	nvec = geomfunc(xy(Ifac,:)); nvec = -nvec(:,2:4);
%    	else
    	nvec = cross(xy(circles(:,2),:)-xy(badnds_fac,:),xy(circles(:,1),:)-xy(badnds_fac,:),2); %points out
    	badbndfac = sort([badnds_fac circles(:,1:2)],2);
    	circID = bndmesh.IDs(bks_bndedg2edg(bndmesh.fac,nd2bndfac,badbndfac));
    	end;
    end;
    if size(xy,2) == 2 || any(Ifac)
    	circles = circles(:,1:max(sum(circles~=0,2)));
    end;
    if size(xy,2) == 2
     if isfield(bndmesh,'triID')
      if any(and(badnds_,bndnds2_))
    	 [circles,badnds,badIDs] = split_2D_circles(circles,badnds,allIDs,bndnds2_,bndmesh,nd2bndedg,nd2tri,tri);
      else
    	badIDs = allIDs(badnds);
      end;
     else
    	badIDs = badnds;
     end;
     [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDt_] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg_);
    else %3D
     if any(Ifac)
      if any(bndnds_(badnds))
      [spheres,newfacs_,newIDs_,newtriN2_,badnds,badbadnds_] = close_spheres(circles,circID,nvec,spheres,badnds,xy,Nmetric,badnds_fac,bndnds2_,bndmesh,nd2bndfac,options);  
      newfacs  = [newfacs ; newfacs_];
      newIDs   = [newIDs  ; newIDs_];
      newtriN2 = [newtriN2; newtriN2_];
      ndelnds(badbadnds_) = true;
      badnds_(badbadnds_) = false;
      end;
     end; %any fac
     if numel(badnds) == 0 %nothing to do here
        continue;
     end;
     %mesh spheres
     if isfield(bndmesh,'triID')
    	badnds_fac2 = find(and(badnds_,bndndsI_));
    	if numel(badnds_fac2) ~= 0
    	[spheres,newfacs_,newIDs_,newtriN2_,badnds,badbadnds_] = split_spheres(spheres,tri,badnds,xy,Nmetric,badnds_fac2,bndndsI_,bndnds2_,bndmesh,nd2bndfac,nd2tri(badnds_fac2,:)',options);
   	newfacs  = [newfacs ; newfacs_];
        newIDs   = [newIDs  ; newIDs_];
    	newtriN2 = [newtriN2; newtriN2_];
    	ndelnds(badbadnds_) = true;
    	end;
    	badIDs = get_badIDs3D(spheres(:,1,:),tri,badnds,bndmesh,nd2tri(badnds,:)'); 
     else
    	badIDs = badnds;
     end;
     [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDt_] = fill_spheres(spheres,Nmetric,xy,badnds,badIDs,triQtb,options);
    end; %triangulation
    
    if numel(badnds) == 0 %we could not improve the situation
        continue;
    end;    
    
    if options.debug == 2 && options.area ~= 0
     %sanity_check3(badnds,newtri_,newtriN_)
     badnds = unique(badnds);
     nd2tri_ = nd2tri(badnds,:);
     [tmp,volo_] = elem_inv(tri(nd2tri_(nd2tri_~=0),:),xy);
     volo = zeros(size(nd2tri_)); volo(nd2tri_~=0) = volo_;
     [newtriN_,I] = sort(newtriN_); newtri_ = newtri_(I,:); qualityN = qualityN(I); newIDt_ = newIDt_(I);
     [tmp,voln_] = elem_inv(newtri_,xy);
     glb2lcl = zeros(max(newtriN_),1);
     if numel(newtriN_) ~= 1
     glb2lcl(newtriN_([newtriN_(1:end-1)~=newtriN_(2:end); true])) = 1:numel(badnds);
     nd2tri_ = inv_table(glb2lcl(newtriN_))';
     else
    	glb2lcl(newtriN_) = 1; nd2tri_ = 1;
     end;
     voln = zeros(size(nd2tri_)); voln(nd2tri_~=0) = voln_;
     if any(abs(sum(voln,1)'-sum(volo,2)) > max(options.minA,1e-12))
    	save for_debug3D.mat; error('domain volume reduced in coarsening, set options.area=0 for curved geometries');
     end;
    end;
    
    %update tri, coordinate and Nmetric table
    ndelnds(badbadnds_) = true; %for active set
    delnds(badnds) = true;
    % fix tri table (new node numbers)
    newtri = [newtri; newtri_];
    newQ = [newQ; qualityN];       
    newIDt = [newIDt; newIDt_];
       
    if i < clrs %update without running gen_books
        %uncolor neighbours to badnds
        clrlssnds = ngh(badnds,:);
        nd2clr(clrlssnds(clrlssnds~=0)) = 0; 
    end;
end;
if any(delnds)
ndone = nnz(delnds);
deltri_ = nd2tri(delnds,:);
I = true(size(tri,1),1); I(deltri_(deltri_~=0)) = false;
tri = [tri(I,:); newtri];
triQ = [triQ(I,:); newQ];
if isfield(bndmesh,'triID')
bndmesh.triID = [bndmesh.triID(I); newIDt];
end;

I = true(size(xy,1),1); I(delnds) = false; xy = xy(I,:); Nmetric = Nmetric(I,:);
if exist('Tmetric','var')
Nmetric = [Nmetric Tmetric(I,:)];
end;
old2new = mvndtri(delnds); tri = old2new(tri); 
if options.mntn_bks
bks.rmnd(ndelnds) = false;
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
bks.mvnd = bks.mvnd(I);
bks.rmnd = bks.rmnd(I);
end;
if size(xy,2) == 3
nI = not(I); nI = nI(newtriN2); 
delbndfac = nd2bndfac(delnds,:); 
I = true(size(bndmesh.fac,1),1); I(delbndfac(delbndfac~=0)) = false;
bndmesh.fac = [bndmesh.fac(I,:); sort(newfacs(nI,:),2)];
bndmesh.IDs = [bndmesh.IDs(I,:); newIDs(nI)];
bndmesh.fac = old2new(bndmesh.fac);
else
bndmesh = update_bndmesh(bndmesh,delnds);
bndmesh.edg = old2new(bndmesh.edg);
end;
bndmesh.crnds = old2new(bndmesh.crnds);

ndone = ndone/size(ngh,1);

if options.debug
   sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;

function [spheres,newfacs,newIDs,newtriN,badnds,badbadnds] = close_spheres(circles,circID,nvec,spheres,badnds,xy,Nmetric,badnds_fac,bndnds2_,bndmesh,nd2bndfac,options)
%fix bad edge nodes
IfacEDG = bndnds2_(badnds_fac);
if any(IfacEDG)
 badnds_edg = badnds_fac(IfacEDG);
[circles2,circID2,nvec2] = split_3D_circles(badnds_edg,circles(bndnds2_(badnds_fac),:),bndmesh,xy,nd2bndfac);
 I1 = 1:sum(IfacEDG); I2 = sum(IfacEDG)+1:size(circles2,1);
 circles(IfacEDG,:) = circles2(I1,:);
 circles = [circles; circles2(I2,:)];
circID(IfacEDG) = circID2(I1);
circID = [circID; circID2(I2)];
nvec(IfacEDG,:) = nvec2(I1,:);
nvec = [nvec; nvec2(I2,:)];

badnds_fac_real = [badnds_fac; badnds_edg];
else
 badnds_fac_real = badnds_fac;
end; %any edge

%%mesh boundaries (close spheres)
nvec = nvec./repmat(sqrt(sum(nvec.^2,2)),1,3);
[newfacs,qualityN,badnds_,newtriN,badbadnds,newRs] = fill_circles(circles,Nmetric,xy,badnds_fac_real,(1:numel(circID))',[],options,nvec);
newIDs = circID(newRs);
if numel(badbadnds) ~= 0
 Ikeep = true(max(badnds),1); Ikeep(badbadnds) = false;
 spheres = spheres(:,:,Ikeep(badnds));
 badnds = badnds(Ikeep(badnds));
end;

%CHECK INVERSION AND AREA?
if options.debug == 2 && numel(newfacs) ~= 0
 nvec = nvec(newRs,:);
 Iedg = [IfacEDG; true(sum(IfacEDG),1)];
 if numel(badbadnds) == 0
 Ikeep = true(max(badnds_fac_real),1);
 end;
I = Ikeep(badnds_fac_real);
sanity_check2(xy,circles(I,:),badnds_fac_real(I),newtriN,nvec,newfacs,options);
end; %debug
if numel(newfacs) ~= 0
[newtriN,I] = sort(newtriN); newIDs = newIDs(I); newfacs = newfacs(I,:);
if numel(newtriN) > 1
d = [newtriN(1:end-1)~=newtriN(2:end); true];
else
d = true;
end;
nN = diff([0; find(d)]); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
spheres2 = zeros(3,nMax,numel(badnds));
glb2lcl = zeros(max(badnds),1); glb2lcl(badnds) = 1:numel(badnds);
spheres2(:,reshape(C(I),size(newfacs,1),1)+(glb2lcl(newtriN)-1)*size(spheres2,2)) = newfacs';
spheres = [spheres spheres2];
end;


function tri = mvndtri(badnds) %reduces node numbers
tri = zeros(size(badnds)); tri(not(badnds)) = 1:nnz(not(badnds));

function badedg = find_badedg(xy,Nmetric,edg,options)
edgnowl = elem_qual(edg,xy,Nmetric,options);
badedg = find(options.Llow > edgnowl);

function rbadnds = badedg2badnds(badedg,edg,nds)
rbadnds = [];
% BAD EDGES TO BAD NODES
if numel(badedg) == 0
    return;
end;
for i = 7+23*(nds(2)==3):-1:1
    badnds = edg(badedg,:); badnds = sort(badnds(:)); badnds_ = badnds(find([badnds(1:end-1)~=badnds(2:end); true]));
    nrbadnds = badnds_(find(diff([0; find([badnds(1:end-1)~=badnds(2:end); true])]) == i)); %we take the nodes connected to 5 bad edges first, then 4 ... 
    %I_ = false(nds,1); I_(nrbadnds) = true; nrbadnds = find(I_); %nrbadnds = nrbadnds(find(not(ismember(nrbadnds,crnds)))); 
    rbadnds = [rbadnds; nrbadnds];
    I_ = false(nds(1),1); I_(rbadnds) = true; edg_ = edg(badedg,:); I = reshape(I_(edg_),size(edg_)); %I = ismember(edg(badedg,:),nrbadnds); 
    badedg = badedg(not(or(I(:,1),I(:,2))));
    if numel(badedg) == 0
        break;
    end;
end;
    
function rbadnds = find_badnds(xy,Nmetric,ngh,options)
n1 = repmat((1:size(xy,1))',1,size(ngh,2));
edg = [n1(ngh~=0) ngh(ngh~=0)];
edgnowl = elem_qual(edg,xy,Nmetric,options);
rbadnds = edg(options.Llow > edgnowl,1);
rbadnds_ = false(size(xy,1),1); 
rbadnds_(rbadnds) = true; %rbadnds_(crnds) = false;
rbadnds = find(rbadnds_);


function bndmesh = update_bndmesh(bndmesh,badnds_)
%badnds_ = false(max([bndmesh.edg(:); badnds']),1); badnds_(badnds) = true;
affedg = badnds_(bndmesh.edg);
affedg_ = or(affedg(:,1),affedg(:,2));
if any(affedg_)
    newedg = bndmesh.edg(not(affedg_),:);
    newIDs = bndmesh.IDs(not(affedg_));
    tmp = bndmesh.edg'; affedg_badnd = tmp(affedg'); affIDs = bndmesh.IDs(affedg_);
    [affedg_badnd,I] = sort(affedg_badnd);
    oldedg = bndmesh.edg(affedg_,:); oldedg = oldedg(I,:); affIDs = affIDs(I);
    I = repmat(affedg_badnd,1,2) ~= oldedg;
    oldedg = oldedg'; 
    newedg = [newedg; sort(reshape(oldedg(I'),2,sum(affedg_)/2)',2)];
    newIDs = [newIDs; affIDs(2:2:end)];
    bndmesh.edg = newedg;
    bndmesh.IDs = newIDs;
end;


function [circles,badnds,badIDs] = split_2D_circles(circles,badnds,allIDs,bndnds2_,bndmesh,nd2bndedg,nd2tri,tri)
badnds2 = badnds(bndnds2_(badnds));
nbad = numel(badnds2);
circles_ = circles(bndnds2_(badnds),:);
ngh1 = bndmesh.edg(nd2bndedg(badnds2,1),:)';
ngh2 = bndmesh.edg(nd2bndedg(badnds2,2),:)';
ngh1 = ngh1(ngh1~=repmat(badnds2',2,1));
ngh2 = ngh2(ngh2~=repmat(badnds2',2,1));
[shft1,R] = find(circles_' == repmat(ngh1',size(circles,2),1));
[shft2,R] = find(circles_' == repmat(ngh2',size(circles,2),1));
NN = sum(circles_~=0,2);
if any(shft1~=1) %shft
[C,R] = find(circles_'); C_ = C-shft1(R)+1;  C1m = C_ < 1;
shft2 = shft2-shft1+1; C2m = shft2 < 1;
C_(C1m) = NN(R(C1m)) + C_(C1m);
shft2(C2m) = NN(C2m) + shft2(C2m);
circles_(R+(C_-1)*nbad) = circles_(R+(C-1)*nbad);
else
[C_,R] = find(circles_');
end;
ocircles = zeros(2*nbad,size(circles,2));
I1 = C_ <= shft2(R);
I2 = C_ >= shft2(R);
ocircles(R(I1)+(C_(I1)-1)*2*nbad) = circles_(R(I1)+(C_(I1)-1)*nbad);
ocircles(R(I2)+(C_(I2)-shft2(R(I2)))*2*nbad+nbad) = circles_(R(I2)+(C_(I2)-1)*nbad);
Rlast = (nbad+1:2*nbad)';
ocircles(Rlast+(NN-shft2+1)*2*nbad) = circles_(:,1);
%determine IDs
triC = sort([repmat(badnds2,2,1) ocircles(:,1:2)],2);
%tri1 = zeros(nbad,size(nd2tri,2)); tri2 = tri1; tri3 = tri1;
%tri1(nd2tri_~=0) = tri(nd2tri_(ndtri_~=0),1);
%tri2(nd2tri_~=0) = tri(nd2tri_(ndtri_~=0),2);
%tri3(nd2tri_~=0) = tri(nd2tri_(ndtri_~=0),3);
nd2tri_ = nd2tri(badnds2,:)';
tri_ = zeros(3,size(nd2tri,2),nbad);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
tri1 = reshape(squeeze(tri_(1,:,:)),size(nd2tri,2),nbad); %reshape active for nbad==1
tri2 = reshape(squeeze(tri_(2,:,:)),size(nd2tri,2),nbad);
tri3 = reshape(squeeze(tri_(3,:,:)),size(nd2tri,2),nbad);
[C1,R] = find(and(and(tri1 == repmat(triC(1:nbad,1)',size(nd2tri,2),1),tri2 == repmat(triC(1:nbad,2)',size(nd2tri,2),1)),tri3 == repmat(triC(1:nbad,3)',size(nd2tri,2),1)));
[C2,R] = find(and(and(tri1 == repmat(triC(nbad+1:end,1)',size(nd2tri,2),1),tri2 == repmat(triC(nbad+1:end,2)',size(nd2tri,2),1)),tri3 == repmat(triC(nbad+1:end,3)',size(nd2tri,2),1)));
triF = [nd2tri_(C1+(R-1)*size(nd2tri,2)); nd2tri_(C2+(R-1)*size(nd2tri,2))];
newIDs = bndmesh.triID(triF);
%write
circles(bndnds2_(badnds),:) = ocircles(1:nbad,:);
circles = [circles; ocircles(nbad+1:end,:)];
badIDs = allIDs(badnds);
badIDs(bndnds2_(badnds)) = newIDs(1:nbad);
badIDs = [badIDs; newIDs(nbad+1:end)];
badnds = [badnds; badnds2];


function [ocircles,circID,nvec] = split_3D_circles(badnds,circles,bndmesh,xy,nd2bndfac)
nbad = numel(badnds);
NN = sum(circles~=0,2);
badnds_ = repmat(badnds,1,size(circles,2));
Ic = circles'~=0;
[Cc,R] = find(Ic);
thmp2 = [2:size(circles,2)+1]'; thmp1 = [0:size(circles,2)-1]'; 
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
fac = zeros(numel(R),3);
fac(:,1) = circles(R+(Cc-1)*nbad);
fac(:,2) = circles(R+(Cr-1)*nbad);
fac(:,3) = badnds(R);
fac = sort(fac,2);
IDs = zeros(size(circles));
IDs(R+(Cc-1)*nbad) = bndmesh.IDs(bks_bndedg2edg(bndmesh.fac,nd2bndfac,fac));
d = and([false(nbad,1) IDs(:,1:end-1)~=IDs(:,2:end)],circles~=0); I = sum(d,2) == 2;
NN1 = zeros(nbad,1);
[shft,R_] = find(d(not(I),:)');
NN1(not(I)) = shft-1;
if any(I)
[shft,R_] = find(d(I,:)');
NN1(I) = shft(2:2:end)-shft(1:2:end);
I = repmat(I',size(circles,2),1); Ic_ = Ic(:,I(1,:)); I = I(Ic);
RI = R(I); CcI = Cc(I); 
shft1 = repmat(shft(1:2:end)'-1,size(circles,2),1); shft1 = shft1(Ic_);
CsI = CcI-shft1; Ip = CsI<1; CsI(Ip) = NN(RI(Ip))+CsI(Ip);
circles(RI+(CsI-1)*nbad) = circles(RI+(CcI-1)*nbad);
IDs(RI+(CsI-1)*nbad) = IDs(RI+(CcI-1)*nbad);
end;
circID = [IDs(:,1); IDs((1:nbad)'+NN1*nbad)];
%circID = [IDs((1:nbad)'+(NN1-1)*nbad); IDs(:,1)];
%nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              %xy(circles((1:nbad)'+(NN1+1)*nbad),:)-xy(badnds,:)], ...
             %[xy(circles(:,1),:)-xy(badnds,:); ... 
              %xy(circles((1:nbad)'+NN1*nbad),:)-xy(badnds,:)],2);
%nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              %xy(circles((1:nbad)'+NN1*nbad),:)-xy(badnds,:)], ...
             %[xy(circles(:,1),:)-xy(badnds,:); ... 
              %xy(circles((1:nbad)'+(NN1-1)*nbad),:)-xy(badnds,:)],2);
nvec =  cross([xy(circles(:,2),:)-xy(badnds,:); ...
              xy(circles(:,1),:)-xy(badnds,:)], ...
             [xy(circles(:,1),:)-xy(badnds,:); ... 
              xy(circles((1:nbad)'+(NN-1)*nbad),:)-xy(badnds,:)],2);
I2 = NN1(R)<Cc; I1 = Cc<=1+NN1(R); 
ocircles = zeros(nbad*2,size(circles,2));
ocircles(R(I1)+(Cc(I1)-1)*nbad*2) = circles(R(I1)+(Cc(I1)-1)*nbad);
ocircles(R(I2)+(Cl(I2)-1)*nbad*2+nbad) = circles(R(I2)+(Cc(I2)-1)*nbad);
ocircles(nbad+1:end,end) = circles(:,1);
%ocircles(nbad+1:end,:) = fix_circles(ocircles(nbad+1:end,:));



function badIDs = get_badIDs3D(spheres_,tri,badnds,bndmesh,nd2tri_)
nbad = numel(badnds);
triC = sort([badnds'; squeeze(spheres_)])';
tri_ = zeros([4 size(nd2tri_,1) nbad]);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
tri1 = reshape(squeeze(tri_(1,:,:)),size(nd2tri_,1),nbad); %reshape active for nbad==1
tri2 = reshape(squeeze(tri_(2,:,:)),size(nd2tri_,1),nbad);
tri3 = reshape(squeeze(tri_(3,:,:)),size(nd2tri_,1),nbad);
tri4 = reshape(squeeze(tri_(4,:,:)),size(nd2tri_,1),nbad);
[C,R] = find(and(and(tri1 == repmat(triC(:,1)',size(nd2tri_,1),1),tri2 == repmat(triC(:,2)',size(nd2tri_,1),1)),and(tri3 == repmat(triC(:,3)',size(nd2tri_,1),1),tri4 == repmat(triC(:,4)',size(nd2tri_,1),1))));
triF = nd2tri_(C+(R-1)*size(nd2tri_,1));
badIDs = bndmesh.triID(triF);

function [spheres_,newfacs_,newIDs_,newtriN2_,badnds,badbadnds] = split_spheres(spheres,tri,badnds,xy,Nmetric,badnds_fac2,bndndsI_,bndnds2_,bndmesh,nd2bndfac,nd2tri_,options)
%save for_debug3D.mat;
nbad = numel(badnds_fac2);
badnds_fac2_ = bndndsI_(badnds);
spheres_ = spheres(:,:,badnds_fac2_);
NN = squeeze(sum(spheres_(1,:,:)~=0,2));
[C,R] = find(squeeze(spheres_(1,:,:)));
spheres4 = zeros(4,numel(R));
if nbad ~= 1
spheres4(1:3,:) = spheres_(1:3,C+(R-1)*size(spheres,2));
else
spheres4(1:3,:) = spheres_(1:3,C+(R-1));
R = C; C=(1:numel(R))';
end;
spheres4(4,:) = badnds_fac2(R);

triC = sort(spheres4);
tri_ = zeros([4 size(nd2tri_,1) nbad]);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
tri1 = reshape(squeeze(tri_(1,:,R)),size(nd2tri_,1),numel(R)); %reshape active for nbad==1
tri2 = reshape(squeeze(tri_(2,:,R)),size(nd2tri_,1),numel(R));
tri3 = reshape(squeeze(tri_(3,:,R)),size(nd2tri_,1),numel(R));
tri4 = reshape(squeeze(tri_(4,:,R)),size(nd2tri_,1),numel(R));
[C_,R_] = find(and(and(tri1 == repmat(triC(1,:),size(nd2tri_,1),1),tri2 == repmat(triC(2,:),size(nd2tri_,1),1)),and(tri3 == repmat(triC(3,:),size(nd2tri_,1),1),tri4 == repmat(triC(4,:),size(nd2tri_,1),1))));
R = reshape(R(R_),numel(R),1); C = C(R_); CR = C+(R-1)*size(spheres,2);
triF = nd2tri_(C_+(R-1)*size(nd2tri_,1));
sphrID = nan(size(spheres_,2),nbad);
sphrID(CR) = bndmesh.triID(triF);
[sphrID,I] = sort(sphrID);
nNaN = not(isnan(sphrID));
spheres_(:,nNaN) = spheres_(:,I(CR)+(R-1)*size(spheres,2));
[C,R] = find(and([false(1,nbad); sphrID(1:end-1,:)~=sphrID(2:end,:)],nNaN));
spheres2 = zeros(3,size(spheres_,2),nbad+numel(R));
badnds2 = zeros(nbad+numel(R),1);
if nbad ~= 1
CR = [ones(nbad,1) rpval2M(R,C) zeros(nbad,1)];
else
CR = [ones(nbad,1) C zeros(nbad,1)];
end;
NN2 = sum(CR~=0,2);
CR((1:nbad)'+NN2*nbad) = NN+1;
ii=1; maxN = max(NN);
Call = repmat((1:maxN)',1,nbad);
Rall = repmat((1:nbad),maxN,1);
for i=1:size(CR,2)-1
   R = find(CR(:,i+1)~=0);
   badnds2(ii:ii-1+numel(R)) = badnds_fac2(R);
   C1 = repmat(CR(R,i  )',maxN,1);
   C2 = repmat(CR(R,i+1)',maxN,1);
   Call_ = Call(:,R); Rall_ = Rall(:,R);
   I = and(C1 <= Call_,Call_ < C2);
   C_ = Call_(I); R_ = Rall_(I);
   Cnew = C_ + 1 - CR(R_,i); 
   if numel(R) == nbad
   Rnew = R_+ii-1;
   else
   R2new = zeros(max(R_),1); tmp = R_([true; R_(1:end-1)~=R_(2:end)]);
   R2new(tmp) = 1:numel(tmp); Rnew = R2new(R_)+ii-1;
   end;
   spheres2(:,Cnew+(Rnew-1)*size(spheres,2)) = spheres_(:,C_+(R_-1)*size(spheres,2));
   ii = ii + numel(R);
end;
I = squeeze(spheres2(1,:,:)) ~= 0;
[C,R] = find(I);
edga = reshape([spheres2(1,I); spheres2(2,I); ...
	        spheres2(2,I); spheres2(3,I); ...
                spheres2(3,I); spheres2(1,I)],2,numel(R)*3);
edga2R   = reshape(repmat(R',3,1),numel(R)*3,1);
[edga,Ia] = sort(edga);
[edga,I] = sortrows([edga2R edga']); 
d = [any(edga(1,:) ~= edga(2,:));
  and(or(edga(2:end-1,1) ~= edga(3:end,1),...
      or(edga(2:end-1,2) ~= edga(3:end,2),...
         edga(2:end-1,3) ~= edga(3:end,3))),...
      or(edga(2:end-1,1) ~= edga(1:end-2,1),...
      or(edga(2:end-1,2) ~= edga(1:end-2,2),...
	 edga(2:end-1,3) ~= edga(1:end-2,3))));...
	 any(edga(end,:) ~= edga(end-1,:))];
edg = edga(d,:);
Iflp = Ia(1,I(d)) == 2;
edg(Iflp,2:3) = edg(Iflp,[3 2]);
Is = find([edg(1:end-1,1)~=edg(2:end,1); true]);
edgN = diff([0; Is]);
circles = zeros(numel(edgN),max(edgN));
circles(:,[1 2]) = edg(Is,[2 3]); %right thumb points in
for i=3:size(circles,2)
   Re = edgN >= i; 
   Re_ = zeros(numel(edgN),1); Re_(Re) = 1:nnz(Re);
   edg_ = edg(Re(edg(:,1)),:);
   ndsrch = circles(Re,i-1);
   Rf = edg_(:,2) == ndsrch(Re_(edg_(:,1)));
   circles(Re,i) = edg_(Rf,3);
end;
nvec = cross(xy(circles(:,2),:)-xy(badnds2,:),xy(circles(:,1),:)-xy(badnds2,:),2);    %points out

circID = bndmesh.IDs(bks_bndedg2edg(bndmesh.fac,nd2bndfac,sort([badnds2 circles(:,1:2)],2)));
Iedg = bndnds2_(badnds2);
badnds_edg = badnds2(Iedg);
[circles2,circID2,nvec2] = split_3D_circles(badnds_edg,circles(Iedg,:),bndmesh,xy,nd2bndfac);
sphr2R = [(1:numel(badnds2))'; find(Iedg)];
badnds2_ = [badnds2; badnds_edg];
circles(Iedg,:) = circles2(1:numel(circID2)/2,:);
nvec(   Iedg,:) = nvec2(   1:numel(circID2)/2,:);
circID(Iedg)    = circID2( 1:numel(circID2)/2);
circles = [circles; circles2(numel(circID2)/2+1:end,:)];
nvec    = [nvec   ; nvec2(   numel(circID2)/2+1:end,:)];
circID  = [circID; circID2(  numel(circID2)/2+1:end)];
save for_debug3D.mat;
[tmp,circR] = sortrows([badnds2_ circID]);
badnds2_ = tmp(:,1); circID = tmp(:,2); circles = circles(circR,:); nvec = nvec(circR,:); sphr2R = sphr2R(circR); circR = (1:numel(circR))';
d2a = [false; and(tmp(1:end-1,1) == tmp(2:end,1),tmp(1:end-1,2) == tmp(2:end,2))];
d2b = [d2a(2:end); false]; nd2 = not(or(d2a,d2b));
d = or(d2a,nd2);
nvec = nvec(d,:)./repmat(sqrt(sum(nvec(d,:).^2,2)),1,3);

[newfacs_,qualityN,badnds_,newtriN2_,badbadnds,new2R] = fill_circles(circles(d,:),Nmetric,xy,badnds2_(d),circR(d),[],options,nvec);
newIDs_ = circID(new2R);
d2a_ = d2a(new2R);
newIDs  = [newIDs_ ; newIDs_(d2a_)];
newfacs = [newfacs_; newfacs_(d2a_,[2 1 3])];
new2R   = sphr2R([new2R; new2R(d2a_)-1]);

if numel(badbadnds) ~= 0
 Ikeep = true(max(badnds),1); Ikeep(badbadnds) = false;
 spheres2 = spheres2(:,:,Ikeep(badnds2));
 badnds2 = badnds2(Ikeep(badnds2));
 spheres = spheres(:,:,Ikeep(badnds));
 badnds  = badnds(Ikeep(badnds));
 old2new = zeros(max(new2R),1); old2new(new2R) = 1; old2new(old2new==1) = 1:sum(old2new);
 new2R = old2new(new2R);
end;


if numel(new2R) ~= 0
   [new2R,I] = sort(new2R); newfacs = newfacs(I,:); newIDs = newIDs(I);
end;
if numel(new2R) > 1
d = [new2R(1:end-1)~=new2R(2:end); true];
else
d = true;
end;
nN = diff([0; find(d)]); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
spheres3 = zeros(3,nMax,numel(badnds2));
spheres3(:,reshape(C(I),size(newfacs,1),1)+(new2R-1)*size(spheres3,2)) = newfacs';
spheres2 = [spheres2 spheres3];

spheres = spheres(:,:,not(bndndsI_(badnds))); badnds = [badnds(not(bndndsI_(badnds))); badnds2];
spheres_ = zeros(3,max(size(spheres,2),size(spheres2,2)),numel(badnds));
spheres_(:,1:size(spheres,2),1:size(spheres,3)) = spheres;
spheres_(:,:,size(spheres,3)+1:end) = [spheres2 zeros(3,size(spheres_,2)-size(spheres2,2),size(spheres2,3))];
%save for_debug3D.mat;

function [] = sanity_check2(xy,circles,badnds,newtriN,nvec,newfacs,options)
circles = fix_circles(circles);
[Ibad,z] = elem_inv(newfacs,xy,nvec); %Ibad = z<options.minA;
if any(Ibad)
	error(sprintf('Coarsening inverted a face trying to close a sphere (%1.1e,%1.1e,%1.1e)',mean(xy(newfacs(find(Ibad,1),:),:))));
end;
nc = size(circles,1); NN = sum(circles~=0,2);
[R,Cc] = find(circles~=0); 
R = reshape(R,numel(R),1); Cc = reshape(Cc,numel(R),1);%reshape when nc == 1
thmp1 = [0:size(circles,2)-1]'; 
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
fac = [reshape(circles(R+(Cc-1)*nc),size(R)) reshape(circles(R+(Cl-1)*nc),size(R)) badnds(R)];
[tmp,zold] = elem_inv(fac,xy);

zoldA = zeros(size(circles)); zoldA(circles~=0) = zold;
a2b = zeros(max(badnds),1); a2b(badnds) = 1:numel(badnds);
newtriN = a2b(newtriN); [newtriN,I] = sort(newtriN); newfacs = newfacs(I,:); z = z(I);
if options.area~=0 && abs(sum(zold)-sum(z)) > max(options.minA,1e-11)
	error(sprintf('domain surface reduced (%1.1e), (curved geometry?)',sum(zold)-sum(z)))
end;

function [circles] = fix_circles(circles)
ocircles = circles';
[Cg,Rg] = find(ocircles~=0);
[C ,R ] = find(ocircles==0);
move = zeros(size(circles)); move(R+(C-1)*size(circles,1)) = -1; move = cumsum(move')';
Cg = Cg + reshape(move(Rg+(Cg-1)*size(circles,1)),size(Cg));
circles = zeros(size(circles));
circles(Rg+(Cg-1)*size(circles,1)) = ocircles(ocircles~=0);
