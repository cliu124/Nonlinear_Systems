function [tri,triQ,bks,bndmesh,ndone] = adapt_flipedg(tri,xy,Nmetric,bndmesh,triQ,bks,options)
% GENERATE BOOKS
ndone = 0;
if size(xy,2) == 2 || options.swap3D
if size(xy,2) == 3
[fac,fac2tri,tri2fac] = bks_all3D(tri);
edg = fac;
nd2edg = inv_table(fac);
edg2tri = fac2tri;
tri2edg = tri2fac;
nbndedg = bks_bndedg2edg(fac,nd2edg,bndmesh.fac);
if size(Nmetric,2) ~= 6
Nmetric = Nmetric(:,1:6);
end;
else
	[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
	nbndedg = bks_bndedg2edg(edg,nd2edg,bndmesh.edg);
	if size(Nmetric,2) ~= 3
	Nmetric = Nmetric(:,1:3);
	end;
end;
%calculate ngh
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
%boundary edges

%color
edg2clr = bks_clr(edg_ngh,options); edg2clr(nbndedg) = 0;
%edgalt = gen_edgalt(edg,tri,edg2tri);
edgalt = gen_edgalt(edg,tri,edg2tri,nbndedg);
clrs = max(edg2clr); %size(edgclr,1);
edgs = size(edgalt,1);
%calculate boundary edge numbers (could be avoided by updating edg instead of regenerating it)

newtri = [];
newQ = [];
newIDs = [];
flipedg = false(size(edg,1),1);

for i=1:clrs
    % FIND EDGES TO FLIP    
    edgi = edg2clr==i;
    if not(any(edgi)) %nothing to do here(should almost never trigger)
    	continue;
    end;
    [edgi,newtri_,Nquality] = find_flipedg(xy,edg,edg2tri,edgalt,edgi,Nmetric,triQ,options);
    if any(edgi) 
            flipedg(edgi) = true;
            newtri = [newtri; newtri_];
            newQ = [newQ; Nquality];
            if isfield(bndmesh,'triID')
            	newIDs = [newIDs; repmat(bndmesh.triID(edg2tri(edgi,1)),size(edg,2),1)];
            end;
            if i<clrs %we are uncolouring the affected edges
                clrlssedg = edg_ngh(edgi,:);
                edg2clr(clrlssedg(clrlssedg ~= 0)) = 0;
            end;
            if options.debug
            	 sanity_check(tri,xy,triQ,Nmetric,options);
            end;
    end;
end;
if any(flipedg)
Nedgi = nnz(flipedg);
ndone = Nedgi;
deltri = reshape(edg2tri(flipedg,:),2*Nedgi,1);
tri(deltri,:) = newtri(1:2*Nedgi,:);
triQ(deltri) = newQ(1:2*Nedgi);
if isfield(bndmesh,'triID')
	bndmesh.triID(deltri) = newIDs(1:2*Nedgi);
end;
if size(edg,2) == 3
tri = [tri; newtri(2*Nedgi+1:end,:)];
triQ = [triQ; newQ(2*Nedgi+1:end)];
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID; newIDs(2*Nedgi+1:end)];
end;
end;
if options.mntn_bks
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
%newedg = edgalt(edgi,:); delfac = find(edgi);
%if size(edg,2) == 3
%newfac = [repmat(newedg,3,1) reshape(edg(edgi,:),Nedgi*3,1)];
%end;
end;
ndone = ndone / size(edg,1);
if options.debug
    sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;
end; %swap3D
if size(xy,2) == 3
	[tri,triQ,bks,bndmesh,ndone_] = adapt_flipedg3D(tri,xy,Nmetric,bndmesh,triQ,bks,options);
	ndone = ndone + ndone_;
end;

function [edgi,newtri,Nquality] = find_flipedg(xy,edg,edg2tri,edgalt,edgi,Nmetric,triQ,options);
dim = size(edg,2);
edgs = nnz(edgi);
if dim == 2	
ntri = [[edgalt(edgi,[2 1]) edg(edgi,1)]; [edgalt(edgi,:) edg(edgi,2)]];
else %3D
ntri = [[edgalt(edgi,:) edg(edgi,[2 1])]; [edgalt(edgi,:) edg(edgi,[3 2])];  [edgalt(edgi,:) edg(edgi,[1 3])]];
end;

qualityO = triQ(edg2tri(edgi,:)); 
qualityOmin = min(reshape(qualityO,edgs,2),[],2);
qualityN = elem_qual(ntri,xy,Nmetric,options,1); 
qualityNmin = min(reshape(qualityN,edgs,dim),[],2);

I = qualityNmin > qualityOmin;
edgi(edgi) = I;
I = repmat(I,dim,1);
newtri = ntri(I,:);
Nquality = qualityN(I,:);


function [tri,triQ,bks,bndmesh,ndone] = adapt_flipedg3D(tri,xy,Nmetric,bndmesh,triQ,bks,options)
%if options.debug
	%sanity_check(tri,xy,triQ,Nmetric,options);
%end;
% GENERATE BOOKS
ndone = 0;
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac,nd2edg,nd2tri] = bks_all3D(tri);
nd2tri = edg2tri;
[tmp,edg1,edg2,edg2fac,edg2ID] = geom_crnds(bndmesh); 
bndedg_ = false(size(edg,1),1); 
bndedg22edg = bks_bndedg2edg(edg,nd2edg,edg2);
bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
bndedg_(bndedg22edg) = true;
if not(isfield(bndmesh,'triID'))
	ngh3 = bks_nd2tri2ndO(tri,nd2tri,bndedg_,edg);
else
	truebndedg = false(size(edg,1),1);
	truebndedg(fac2edg(fac2tri(:,2)==0,:)) = true;
	ngh3 = bks_nd2tri2ndO(tri,nd2tri,truebndedg,edg);
end;
ngh = bks_edg2tri2edg(edg2tri,tri2edg);
%3D: Do not flip edges on the boundary of the boundary

nd2clr = bks_clr(ngh,options); nd2clr(bndedg2edg) = 0;
if isfield(options,'TS')
	bndfac2fac = bks_bndedg2edg(fac,nd2fac,bndmesh.fac);
	nd2clr(fac2edg(bndfac2fac(bndmesh.IDs==options.TS),:)) = 0;
end;
%if isfield(bndmesh,'triID')  %filthy hack
	%nd2clr(and(not(truebndedg),bndedg_)) = 0;
%end;

triQtb = inf(size(nd2tri)); I = nd2tri~=0;
triQtb(I) = triQ(nd2tri(I)); triQtb = min(triQtb,[],2);
  
newtri = [];
newIDs = [];
newQ = [];
delnds = false(size(nd2tri,1),1); 
ndelnds = false(size(nd2tri,1),1);
clrs = max(nd2clr);
for i=1:clrs
    badnds_ = nd2clr==i; 
    badnds  = find(badnds_);
    if numel(badnds) == 0 %nothing to do here
        continue;
    end;
    circles = ngh3(badnds,:);
    circles = circles(:,1:max(sum(circles~=0,2)));
    if isfield(bndmesh,'triID')
    	badnds2_ = and(not(truebndedg(badnds_)),bndedg_(badnds_));
    	if any(badnds2_)
    		[circles,badnds] = split_3D_circles(circles,tri,edg(badnds,:),badnds,badnds2_,bndmesh,nd2tri(badnds,:));
    	end;
    	badIDs = get_badIDs(circles,tri,edg(badnds,:),bndmesh,nd2tri(badnds,:)');
    else
    	badIDs = badnds;
    end;
    [newtri_,qualityN,badnds,newtriN_,badbadnds_,newIDs_] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg(badnds,:));
    
    if numel(badnds) == 0 %we could not improve the situation
        continue;
    end;    
    %update tri, coordinate and Nmetric table
    ndelnds(badbadnds_) = true;
    delnds(badnds) = true;
    % fix tri table (new node numbers)
    newIDs = [newIDs; newIDs_];
    newtri = [newtri; newtri_];
    newQ = [newQ; qualityN];       
       
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
if isfield(bndmesh,'triID');
bndmesh.triID = [bndmesh.triID(I); newIDs];
end;

bndmesh = update_bndmesh_3D(bndmesh,delnds,edg2fac,edg2,edg2ID, bndedg2edg,bndedg22edg);

if options.mntn_bks
bks.mvnd(newtri) = true;
bks.rmnd(newtri) = true;
end;
ndone = ndone/size(ngh,1);

if options.debug
   sanity_check(tri,xy,triQ,Nmetric,options);
end;
end;


function bndmesh =  update_bndmesh_3D(bndmesh,deledg_,edg2fac,edg2,edg2ID,bndedg2edg,bndedg22edg)
%deledg_ = false(max(max(deledg),max(bndedg22edg)),1); deledg_(deledg) = true;
flipedg = deledg_(bndedg22edg);
if not(any(flipedg))
	return;
end;
DMcnv = zeros(max(edg2(:)),1); DMcnv(edg2(:)) = 1; DMcnv(DMcnv~=0) = 1:nnz(DMcnv);
I_ = false(size(deledg_,1),1); I_(bndedg2edg) = true;
edgalt = gen_edgalt(DMcnv(edg2),DMcnv(bndmesh.fac),edg2fac,find(I_(bndedg22edg)));
DMcnvINV = find(DMcnv);
edgalt(edgalt(:,1)~=0,:) = DMcnvINV(edgalt(edgalt(:,1)~=0,:));
flipfac = edg2fac(flipedg,[1 2]);
newfac = [[edg2(flipedg,1) edgalt(flipedg,:)]; [edg2(flipedg,2) edgalt(flipedg,:)]];
newIDs = repmat(edg2ID(flipedg,1),2,1);
bndmesh.fac(flipfac(:),:) = sort(newfac,2);
bndmesh.IDs(flipfac(:),:) = newIDs;


function [circles,badnds] = split_3D_circles(circles,tri,edg,badnds,badnds2_,bndmesh,nd2tri_);
badnds2 = badnds(badnds2_);
nd2tri_ = nd2tri_(badnds2_,:)';
nbad = numel(badnds2);
edg = edg(badnds2_,:);
circles_ = circles(badnds2_,:);
NN = sum(circles_~=0,2);
badnds_ = repmat(badnds2,1,size(circles,2));
Ic = circles_'~=0;
[Cc,R] = find(Ic);
thmp2 = [2:size(circles_,2)+1]'; thmp1 = [0:size(circles_,2)-1]'; 
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
triC = zeros(numel(R),4);
triC(:,1) = circles_(R+(Cc-1)*nbad); triC(:,3) = edg(R,1);
triC(:,2) = circles_(R+(Cr-1)*nbad); triC(:,4) = edg(R,2);
triC = sort(triC,2);
tri_ = zeros([4 size(nd2tri_,1) nbad]);
tri_(:,nd2tri_~=0) = tri(nd2tri_(nd2tri_~=0),:)';
tri_ = sort(tri_);
%tri1 = reshape(squeeze(tri_(1,:,R)),size(nd2tri_,1),numel(R)); %reshape active for nbad==1
tri1 = squeeze(tri_(1,:,R)); tri2 = squeeze(tri_(2,:,R)); 
tri3 = squeeze(tri_(3,:,R)); tri4 = squeeze(tri_(4,:,R));
[C,R_] = find(and(and(tri1 == repmat(triC(:,1)',size(nd2tri_,1),1),tri2 == repmat(triC(:,2)',size(nd2tri_,1),1)),and(tri3 == repmat(triC(:,3)',size(nd2tri_,1),1),tri4 == repmat(triC(:,4)',size(nd2tri_,1),1))));
triF = nd2tri_(C+(R(R_)-1)*size(nd2tri_,1));
IDs = zeros(size(circles_)); 
IDs(R+(Cc-1)*nbad) =  bndmesh.triID(triF);
d = and([false(nbad,1) IDs(:,1:end-1)~=IDs(:,2:end)],circles_~=0); I = sum(d,2) == 2;
NN1 = zeros(nbad,1);
[shft,R_] = find(d(not(I),:)');
NN1(not(I)) = shft-1;
if any(I)
[shft,R_] = find(d(I,:)');
NN1(I) = shft(2:2:end)-shft(1:2:end);
I = repmat(I',size(circles_,2),1); Ic_ = Ic(:,I(1,:)); I = I(Ic);
RI = R(I); CcI = Cc(I); 
shft1 = repmat(shft(1:2:end)'-1,size(circles_,2),1); shft1 = shft1(Ic_);
CsI = CcI-shft1; Ip = CsI<1; CsI(Ip) = NN(RI(Ip))+CsI(Ip);
circles_(RI+(CsI-1)*nbad) = circles_(RI+(CcI-1)*nbad);
IDs(RI+(CsI-1)*nbad) = IDs(RI+(CcI-1)*nbad);
end;
I2 = NN1(R)<Cc; I1 = Cc<=1+NN1(R); 
ocircles = zeros(nbad*2,size(circles_,2));
ocircles(R(I1)+(Cc(I1)-1)*nbad*2) = circles_(R(I1)+(Cc(I1)-1)*nbad);
ocircles(R(I2)+(Cl(I2)-1)*nbad*2+nbad) = circles_(R(I2)+(Cc(I2)-1)*nbad);
ocircles(nbad+1:end,end) = circles_(:,1);
ocircles = fix_circles(ocircles);
badnds = [badnds; badnds2];
circles(badnds2_,:) = ocircles(1:nbad,:);
circles = [circles; ocircles(nbad+1:end,:)];

function badIDs = get_badIDs(circles,tri,badedg,bndmesh,nd2tri_)
nbad = size(badedg,1);
triC = sort([badedg circles(:,[1 2])],2);
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
if abs(sum(zold)-sum(z)) > options.minA
	error(sprintf('domain surface reduced (%1.1e), (curved geometry?)',sum(zold)-sum(z)))
end;

function circles = fix_circles(circles)
ocircles = circles';
[Cg,Rg] = find(ocircles~=0);
[C ,R ] = find(ocircles==0);
move = zeros(size(circles)); move(R+(C-1)*size(circles,1)) = -1; move = cumsum(move')';
Cg = Cg + reshape(move(Rg+(Cg-1)*size(circles,1)),size(Cg));
circles = zeros(size(circles));
circles(Rg+(Cg-1)*size(circles,1)) = ocircles(ocircles~=0);