function [xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
%GENERATE BOOKS
if size(xy,2) == 3
	[xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd3D(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options);
	return;
end;

ndone = 0;
[edg,edg2tri,tri2edg,nd2edg,edga,edga2tri,tri2edga,edga2edg] = bks_all(tri);
bndedg2edg = bks_bndedg2edg(edg,nd2edg,bndmesh.edg);
[edgnowL,c] = elem_qual(edg,xy,Nmetric(:,1:3),options);
splitedg = options.Lup<edgnowL;
if not(options.fastRFN)
	splitedg = true(size(splitedg));
end;

if not(any(splitedg)) 
    return
end;
if 0 < options.qualP
[triangle,badedg,badangle] = elem_angle(tri,xy,options);
badedg2edg = bks_bndedg2edg(edg,nd2edg,badedg);
splitedg = false(size(splitedg));
splitedg(badedg2edg) = true;
allangle = zeros(size(splitedg)); allangle(badedg2edg) = badangle;
edgnowL = allangle;
end;
if options.smpRFN == 1%split a single edge per element
   splitedg = splitedg1(splitedg,edgnowL,tri2edg,edg2tri);
end;
if isfield(options,'TS')
	splitedg(bndedg2edg(bndmesh.IDs==options.TS)) = false;
end;

if options.smpRFN == 2 %using colouring
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,xyf,ndone] = calc_new_wrap_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
else
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
xyf = find(splitedg);
end;
Nsplt = nnz(splitedg);
if Nsplt==0
	return;
end;

newnds = size(xy,1)+(1:sum(splitedg))';


if options.mntn_bks
bks.mvnd(tri(deltric~=0,:)) = true;
bks.rmnd(tri(deltric~=0,:)) = true;
bks.mvnd =  [bks.mvnd; repmat(true,Nsplt,1)];
bks.rmnd =  [bks.rmnd; repmat(true,Nsplt,1)];
end;

xy = nxy;
Nmetric = nNmetric;
%newtri = elem_fix(newtri,xy);
tri = [tri(deltric==0,:); newtri];
triQ = [triQ(deltric==0); q3]; 
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID(deltric==0); bndmesh.triID(oldtri)];
end;

% UPDATE BOUNDARY MESH all(bndmesh.edg == sort(bndmesh.edg,2)) == true
%calculate boundary edge numbers (could be avoided by updating edg instead of regenerating it)
edgs = size(edg,1);
% delete split edges and add new ones
splitedgbnd = splitedg(bndedg2edg);
if any(splitedgbnd)
    %newnds = size(xy,1)-sum(splitedg)+(1:sum(splitedg))';
    edg2newnd = zeros(edgs,1); edg2newnd(xyf) = newnds;
    newnds = edg2newnd(bndedg2edg(splitedgbnd));
    newedg = bndmesh.edg(not(splitedgbnd),:);
    newIDs = bndmesh.IDs(not(splitedgbnd ));
    newedg = [newedg; [edg(bndedg2edg(splitedgbnd),1) newnds; edg(bndedg2edg(splitedgbnd),2) newnds]];
    newIDs = [newIDs; repmat(bndmesh.IDs(splitedg(bndedg2edg)),2,1)];
    bndmesh.edg = newedg;
    bndmesh.IDs = newIDs;
end;


if options.debug
    sanity_check(tri,xy,triQ,Nmetric(:,1:3),options);
end;

function [newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,xyf,ndone] = calc_new_wrap_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options)

edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
edg2clr = bks_clr(edg_ngh,options); clrs = max(edg2clr);
ndone = 0; newedg = zeros(0,2); oldtri = []; nNmetric = zeros(0,3); nxy = zeros(0,2); q3 = []; newtri = zeros(0,3); splitedgO = false(size(edg,1),1); deltric = zeros(size(tri,1),1); newnds = []; xyf = []; Nsplt = 0;
for i=1:clrs
	splitedg_ = and(splitedg,edg2clr==i);
	[newtri_,q3_,nxy_,nNmetric_,oldtri_,deltric_,newedg_,splitedg_,ndone_] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg_,c,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
	splitedgO(splitedg_) = true;
	fxrr = [(1:size(xy,1)) Nsplt+size(xy,1)+(1:nnz(splitedg_))]';
	newtri_ = fxrr(newtri_); newedg_ = fxrr(newedg_);
	Nsplt = Nsplt+nnz(splitedg_);
	ndone = ndone + ndone_;
	newtri = [newtri; newtri_];
	q3 = [q3; q3_];
	nxy = [nxy; nxy_(size(xy,1)+1:end,:)];
	xyf = [xyf; find(splitedg_)];
	nNmetric = [nNmetric; nNmetric_(size(xy,1)+1:end,:)];
	newedg = [newedg; newedg_];
	oldtri = [oldtri; oldtri_];
	
	deltric(deltric_~=0) = 1;
	if i ~= clrs
	edg2clr(tri2edg(oldtri_,:)) = 0;
	end;
end; %for
nxy = [xy; nxy];
nNmetric = [Nmetric; nNmetric];
splitedg = splitedgO; 

function [newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,tri2edg,tri2edga,bndedg2edg,geomfunc,options)
while true
[newtri,q3,nxy,nNmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,edg,edga,tri2edg,tri2edga,bndedg2edg,geomfunc,options);
badE = [];
goodE = true(size(q3,1),1);
if options.consRFN == 1
	goodE = q3 >= triQ(oldtri);
elseif options.minA ~= 0 || numel(geomfunc)~=0 || options.debug
	goodE = q3 >= options.minqual;
end;
if options.consRFN == 2 && options.smpRFN
badE = calc_badE(oldtri,q3,triQ,splitedg,tri2edg,edg2tri); 
else
badE = oldtri(not(goodE));
end;
if numel(badE) ~= 0 && not(options.consRFN)
	warning(sprintf('RECOMPUTING REFINEMENT due to inverted element (%1.1e,%1.1e), you might want to help us by inserting fixed nodes on your curved concave boundary.',mean(nxy(newtri(find(q3<0,1),:),:))));
end;
	
if numel(badE) ~= 0
	badedg = tri2edg(badE,:);
	if options.smpRFN
	splitedgo = splitedg;
	splitedg(badedg(:)) = false;
	splitedg_ = splitedg(splitedgo);
	nNmetric = nNmetric([true(size(xy,1),1); splitedg_],:);
	nxy = nxy([true(size(xy,1),1); splitedg_],:);
	newtri = newtri(repmat(splitedg_,2,1),:);
	oldtri = oldtri(repmat(splitedg_,2,1),:);
	q3 = q3(repmat(splitedg_,2,1));
	newedg = newedg(repmat(splitedg_,3,1),:);
	deltric = sum(splitedg(tri2edg),2);
	ndone = sum(splitedg)/size(edg,1);
	else
	splitedg(badedg(:)) = false;
	end;
else
	break;
end;
end;

function badE = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri)
edg2tri = edg2tri(splitedg,:);
edgQ = inf(size(edg2tri)); edgQ(edg2tri~=0) = triQ(edg2tri(edg2tri~=0));
edgQ_ = zeros(size(splitedg)); edgQ_(splitedg) = min(edgQ,[],2);
newq_ = min(reshape(newq,size(newq,1)/2,2),[],2);
oldtri = oldtri(1:size(oldtri,1)/2);
tri2edg = tri2edg(oldtri,:)'; edg = tri2edg(splitedg(tri2edg));
Ibad = newq_ < edgQ_(edg);
badE = oldtri(Ibad);


function [newtri,q3,xy,Nmetric,oldtri,deltric,newedg,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,edg,edga,tri2edg,tri2edga,bndedg2edg,geomfunc,options)

newcoords = c(splitedg,:);

%fix newcoords in case of curved geomtry
if numel(geomfunc)~=0
	bndedg = false(size(edg,1),1);
	bndedg(bndedg2edg) = true;
	I = bndedg(splitedg);
	dist = geomfunc(newcoords(I,:));
	newcoords(I,:) = newcoords(I,:) - dist(:,2:3).*repmat(dist(:,1),1,2);
	% we will not allow interior nodes outside the geometry:
	nI = not(I);
	dist = geomfunc(newcoords(nI,:));
	goodI = true(size(I)); goodI(nI) = dist(:,1) > -options.geomtol;
	if any(not(goodI))
		warning(sprintf('REFINEMENT CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(goodI))));
		splitedg(splitedg) = goodI;
		newcoords = newcoords(goodI,:);
	end;
end;
ndone = sum(splitedg)/size(edg,1);

%% disp(sprintf('### MESH REFINEMENT: %0.1f%% more nodes',size(newcoords,1)));

I_ = find(splitedg); 
I2 = zeros(size(edg,1),1); I2(splitedg) = 1:numel(I_);
deltric = sum(splitedg(tri2edg),2);
deltri = I2(tri2edg);
deltri_ = deltri+size(xy,1);
%deltrit = deltri'; deltrit(deltrit==0) = numel(I_)+1; I__ = [I_;0];

%deltri contains the splitedg numbers for each element
%deltri_ contains the new node numbers for each element

%UPDATE Nmetric and calculate new coordinates
ndsO = size(xy,1);
xy = [xy; newcoords];
if size(Nmetric,2) ~= 3
Tmetric = Nmetric(:,4:end); Nmetric = Nmetric(:,1:3); 
Tmetric = [Tmetric; (Tmetric(edg(splitedg,1),:) + Tmetric(edg(splitedg,1),:))/2.];
end;
Nmetric = [Nmetric; metric_avg(edg(splitedg,:),Nmetric,options)];
	 	
% SPLIT ELEMENTS WITH 3 LONG EDGES
I = deltric==3;
if any(I)
    N3 = nnz(I);
    if not(options.advRFN)
    	newtri = [deltri_(I,:); ...
		[tri(I,2) deltri_(I,[2 1])]; ...
		[tri(I,3) deltri_(I,[3 2])]; ...
	        [tri(I,1) deltri_(I,[1 3])]];
        q3 = elem_qual(newtri,xy,Nmetric,options);
    else   	
    nodes = [deltri_(I,:) tri(I,:)]; 
    chc = [1 2 3 5 2 1 6 3 2 4 1 3; ...
	   5 2 3 5 3 1 6 3 2 4 1 3; ...
    	   6 3 1 5 2 1 6 1 2 4 1 3; ...
    	   4 1 2 5 2 1 6 3 2 4 2 3];
    newtri = zeros(4*N3,3); q3 = zeros(size(newtri,1),1); qmin = zeros(N3,1);
    for i=1:4
	 newtri_ = [nodes(:,chc(i,[1 2 3])); nodes(:,chc(i,[4 5 6])); nodes(:,chc(i,[7 8 9])); nodes(:,chc(i,[10 11 12]))];
	 q3_ = elem_qual(newtri_,xy,Nmetric,options);
	 qmin_ = min(reshape(q3_,N3,4),[],2);
	 better = qmin < qmin_;
	 qmin(better) = qmin_(better);
	 better4 = repmat(better,4,1);
	 newtri(better4,:) = newtri_(better4,:);
	 q3(better4) = q3_(better4);
    end;%fori
    end;
    oldtri = repmat(find(I),4,1); 
    newedg = [[tri(I,2) deltri_(I,1)]; [tri(I,2) deltri_(I,2)]; ...
              [tri(I,3) deltri_(I,2)]; [tri(I,3) deltri_(I,3)]; ...
              [tri(I,1) deltri_(I,1)]; [tri(I,1) deltri_(I,3)]; ...
              newtri([false(N3,1); true(3*N3,1)],[2 3])];
else
    newtri = zeros(0,3); q3 = []; oldtri = []; newedg = zeros(0,2);
end;


% SPLIT ELEMENTS WITH 1 LONG EDGE
I = deltric==1;
if any(I)
    [C1,R_] = find(deltri(I,:)'~=0);
    C1a = C1+1; C1a(C1a==4) = 1; C1b = C1-1; C1b(C1b==0) = 3; R = find(I);
    oppnd = tri(R+size(tri,1)*(C1b-1));
    ooppnd = deltri_(R+size(tri,1)*(C1-1));
    newtri2 = [[oppnd tri(R+size(tri,1)*(C1-1)) ooppnd]; ...
               [oppnd ooppnd tri(R+size(tri,1)*(C1a-1))]];
    newtri = [newtri; newtri2];     
    q3 = [q3; elem_qual(newtri2,xy,Nmetric,options)];
    oldtri = [oldtri; repmat(R,2,1)];
    newedg = [newedg; [ooppnd oppnd]; newtri2(:,[2 3])];
end; 

% SPLITE ELEMENTS WITH 2 LONG EDGES
I = deltric==2;
if any(I)    
    [C2,R] = find(deltri(I,:)'==0);
R = find(I);
    C2b = C2-1; C2b(C2b==0) = 3;
    shrnd = tri(R+(C2b-1)*size(tri,1));
    edga2 = tri2edga(R+(C2-1)*size(tri,1));   
    deltri2_ = deltri_(I,:)';
    nshrnd = reshape(deltri2_(deltri2_~=ndsO),2,size(deltri2_,2))';
    tswitch = C2 ~= 2;
    nshrnd(tswitch,:) = nshrnd(tswitch,[2 1]);
    edgO = [edga(edga2,:) nshrnd];
    nI2 = size(edgO,1);
    tri1 = [edgO(:,[2 1 4]); edgO(:,[2 4 3])];
    tri2 = [edgO(:,[3 1 4]); edgO(:,[1 3 2])]; 
    %% DETERMINE HOW TO SPLIT
    quality1 = elem_qual(tri1,xy,Nmetric,options); 
    quality1m = min(reshape(quality1,nI2,2),[],2);
    quality2 = elem_qual(tri2,xy,Nmetric,options); 
    quality2m = min(reshape(quality2,nI2,2),[],2);
    option1_ = quality2m<quality1m; 
    option1 = repmat(option1_,2,1);
    option2 = not(option1); 
    newtri3 = [[shrnd nshrnd]; tri1(option1,:); tri2(option2,:)];
    newtri = [newtri; newtri3];
    q3 = [q3; elem_qual([shrnd nshrnd],xy,Nmetric,options); quality1(option1); quality2(option2)];
    R = find(I);
    oldtri = [oldtri; R; repmat(R(option1_),2,1); repmat(R(not(option1_)),2,1)];
    newedg = [newedg; nshrnd; newtri3([false(2*nI2,1); true(nI2,1)],[1 2])];
end;
if exist('Tmetric','var')
Nmetric = [Nmetric Tmetric];
end;


function splitedg_ = splitedg1(splitedg,edgnowL,tri2edg,edg2tri)
edgs = size(splitedg,1);
splitedg_ = false(edgs,1);
splitedgf = zeros(edgs,1); splitedgf(splitedg) = find(splitedg);
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
nghI = edg_ngh~=0;
while any(splitedg)
edgnowL(not(splitedg)) = 0;
cmp = zeros(size(edg_ngh)); cmp(nghI) = edgnowL(edg_ngh(nghI));
cmp2 = zeros(size(edg_ngh)); cmp2(nghI) = splitedgf(edg_ngh(nghI));
I = all(or(repmat(edgnowL,1,size(edg_ngh,2)) > cmp,and(and(repmat(edgnowL,1,size(edg_ngh,2))~=0,repmat(edgnowL,1,size(edg_ngh,2)) == cmp),repmat(splitedgf,1,size(edg_ngh,2))>cmp2)),2); %sum(I)
splitedg_(I) = true;
nsplit = edg_ngh(I,:); nsplit = nsplit(nsplit(:)~=0);
splitedg(nsplit(:)) = false; splitedg(I) = false; splitedgf(I) = 0; 
end;