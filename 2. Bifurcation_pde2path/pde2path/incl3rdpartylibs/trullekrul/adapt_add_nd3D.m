function [xy,tri,bndmesh,Nmetric,triQ,bks,ndone] = adapt_add_nd3D(tri,xy,Nmetric,bndmesh,triQ,bks,geomfunc,options)
%GENERATE BOOKS
[fac,fac2tri,tri2fac,faca,faca2tri,tri2faca,faca2fac,edg,edga,edga2edg,edg2tri,tri2edg,tri2edga,edga2tri,fac2edg,edg2fac,edga2faca,faca2edga,faca2edg,fac2edga,fac2tri2,nd2fac] = bks_all3D(tri);
bndfac2fac = bks_bndedg2edg(fac,nd2fac,bndmesh.fac);
ndone = 0;
[edgnowL,c] = elem_qual(edg,xy,Nmetric(:,1:6),options);
splitedg = options.Lup<edgnowL;
if not(options.fastRFN)
	splitedg = true(size(splitedg));
end;
if not(any(splitedg)) 
    return;
end;
if 0 < options.qualP
nd2edg = inv_table(edg);
[triangle,badedg,badangle] = elem_angle(tri,xy,options);
badedg2edg = bks_bndedg2edg(edg,nd2edg,badedg);
splitedg = false(size(splitedg));
splitedg(badedg2edg) = true;
allangle = zeros(size(splitedg)); allangle(badedg2edg) = badangle;
edgnowL = allangle;
end;
if numel(geomfunc) ~= 0
    nd2edg = inv_table(edg);
    [crnds,edg1,edg2] = geom_crnds(bndmesh); bndedg2edg = bks_bndedg2edg(edg,nd2edg,edg1);
else
    bndedg2edg = [];
end;

if options.smpRFN == 1 %split a single edge per element
   splitedg = splitedg1(splitedg,edgnowL,tri2edg,edg2tri);
end;
if isfield(options,'TS')
	splitedg(fac2edg(bndfac2fac(bndmesh.IDs==options.TS),:)) = false;
end;

if options.smpRFN == 2
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg);
edg2clr = bks_clr(edg_ngh,options); clrs = max(edg2clr);
ndone = 0; newfacs = zeros(0,3); newIDs = []; oldtri = []; nNmetric = zeros(0,size(Nmetric,2)); nxy = zeros(0,3); newq = []; newtri = zeros(0,4); splitedgO = false(size(edg,1),1); deltric = zeros(size(tri,1),1); keepbndfacs=true(size(bndmesh.IDs)); Nsplt = 0;
for i=1:clrs
	splitedg_ = and(splitedg,edg2clr==i);
        [newtri_,newq_,nxy_,nNmetric_,oldtri_,deltric_,keepbndfacs_,newfacs_,newIDs_,splitedg_,ndone_] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg_,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
	splitedgO(splitedg_) = true;
	fxrr = [(1:size(xy,1)) Nsplt+size(xy,1)+(1:nnz(splitedg_))]';
	newtri_ = fxrr(newtri_); newfacs_ = fxrr(newfacs_);
	Nsplt = Nsplt+nnz(splitedg_);
	ndone = ndone + ndone_;
	newtri = [newtri; newtri_];
	newq = [newq; newq_];
	nxy = [nxy; nxy_(size(xy,1)+1:end,:)];
	nNmetric = [nNmetric; nNmetric_(size(xy,1)+1:end,:)];
	newfacs = [newfacs; newfacs_];
	newIDs = [newIDs; newIDs_];
	oldtri = [oldtri; oldtri_];
	keepbndfacs(not(keepbndfacs_)) = false;
	deltric(deltric_~=0) = 1;
	if i ~= clrs
	edg2clr(tri2edg(oldtri_,:)) = 0;
	end;
end; %for
nxy = [xy; nxy];
nNmetric = [Nmetric; nNmetric];
splitedg = splitedgO; 
else
[newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
end;
if not(any(splitedg))
	return;
end;
xy = nxy;
Nmetric = nNmetric;
bndmesh.fac = [bndmesh.fac(keepbndfacs,:); sort(newfacs,2)];
bndmesh.IDs = [bndmesh.IDs(keepbndfacs); newIDs];
if isfield(bndmesh,'triID')
	bndmesh.triID = [bndmesh.triID(deltric==0); bndmesh.triID(oldtri)];
end;
%trio=tri;

if options.mntn_bks && any(splitedg)
bks.mvnd(tri(deltric~=0,:)) = true;
bks.rmnd(tri(deltric~=0,:)) = true;
Nsplt = size(xy,1)-size(nd2fac,1);
bks.mvnd =  [bks.mvnd; true(Nsplt,1)];
bks.rmnd =  [bks.rmnd; true(Nsplt,1)];
end;

tri = [tri(deltric==0,:); newtri];
triQ = [triQ(deltric==0); newq]; 

if options.debug
sanity_check(tri,xy,triQ,Nmetric(:,1:6),options);
if numel(geomfunc) ~= 0 && size(xy,2) == 3
bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.fac) = true;
dist = geomfunc(xy(bndnds_,:));
if max(abs(dist(:,1))) > options.geomtol
	error('boundary node far from boundary');
end;
end;
end;

function [newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = calc_new_wrap(tri,xy,Nmetric,triQ,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options)


while true
[newtri,newq,nxy,nNmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options);
goodE = true(size(newq,1),1); 
if options.consRFN == 1
	goodE = newq >= triQ(oldtri);
elseif options.minA ~= 0 || numel(geomfunc)~=0 || options.debug
	goodE = newq >= options.minqual;
end; 
if options.consRFN == 2 && options.smpRFN
badE = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri); 
%[badE,goodE,splitedg,goodnd] = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri); 
%deltric(oldtri(not(goodE))) = 0; oldtri = oldtri(goodE); newq = newq(goodE);  old2new = zeros(size(goodnd)); old2new(goodnd) = 1:sum(goodnd); old2new =  [(1:size(xy,1))'; old2new+size(xy,1)]; newtri = old2new(newtri(goodE,:)); goodnd = [true(size(xy,1),1); goodnd]; nxy = nxy(goodnd,:); nNmetric = nNmetric(goodnd,:); goodfacs = all(goodnd(newfacs),2); newfacs = old2new(newfacs(goodfacs,:)); newIDs = newIDs(goodfacs);  keepbndfacs = true(size(bndfac2fac)); keepbndfacs(any(splitedg(fac2edg(bndfac2fac,:)),2)) = false; ndone = sum(splitedg)/size(edg,1); break;
else
badE = oldtri(not(goodE));
end;

if numel(badE) ~= 0 && not(options.consRFN)
	warning(sprintf('RECOMPUTING REFINEMENT due to inverted element (%1.1e,%1.1e,%1.1e), you might want to help us by inserting fixed nodes on your curved boundary.',mean(nxy(newtri(find(newq<0,1),:),:))));
end;

if numel(badE) ~= 0
	badedg = tri2edg(badE,:);
	if options.smpRFN
		deledg = false(size(splitedg)); deledg(badedg) = true; deledg(not(splitedg)) = false;
		edg2tri_ = edg2tri(deledg,:); ndeltri = true(size(tri,1),1);
		ndeltri(edg2tri_(edg2tri_~=0)) = false;
		newtri = newtri(ndeltri(oldtri),:);
		newq   =   newq(ndeltri(oldtri));
		oldtri = oldtri(ndeltri(oldtri));
		deltric(not(ndeltri)) = 0;
		I = [true(size(xy,1),1); not(deledg(splitedg))];
		nNmetric = nNmetric(I,:);
		old2new = zeros(size(nxy,1),1);
		nxy = nxy(I,:);
		old2new(I) = 1:size(nxy,1);
		newIDs  =  newIDs(I(newfacs(:,1)));
		newfacs = newfacs(I(newfacs(:,1)),:);
		keepbndfacs(any(deledg(fac2edg(bndfac2fac,:)),2)) = true;
		splitedg(deledg) = false;
		newtri(:,1)  = old2new(newtri(:,1));
		newfacs(:,1) = old2new(newfacs(:,1));
		ndone = sum(splitedg)/size(edg2fac,1);
		break;
	else
		splitedg(badedg(:)) = false;
	end;
else
	 break; 
end;
end; %while


function [newtri,newq,xy,Nmetric,oldtri,deltric,keepbndfacs,newfacs,newIDs,splitedg,ndone] = calc_new(tri,xy,Nmetric,splitedg,c,edg,edga,edg2tri,fac,tri2fac,tri2edg,tri2edga,edg2fac,fac2edg,fac2edga,fac2tri2,bndmesh,bndedg2edg,bndfac2fac,geomfunc,options)
newcoords = c(splitedg,:);
%fix newcoords in case of curved geomtry
if numel(geomfunc)~=0 
	bndedg  = false(size(edg,1),1);
 	bndedg2 = false(size(edg,1),1);
	bndedg(fac2edg(bndfac2fac,:)) = true; 
	I = bndedg(splitedg);
	bndedg(bndedg2edg) = false;
	bndedg2(bndedg2edg) = true;
	I1 = bndedg(splitedg);
	I2 = bndedg2(splitedg);
	dist = zeros(sum(I),4);
	dist(I1(I),:) = geomfunc{1}(newcoords(I1,:));
	if numel(bndedg2edg) ~= 0
	dist(I2(I),:) = geomfunc{2}(newcoords(I2,:));
	end;
	newcoords(I,:) = newcoords(I,:) - dist(:,2:4).*repmat(dist(:,1),1,3);	
	% we will not allow interior nodes outside the geometry:
	nI = not(I);
	dist = geomfunc{1}(newcoords(nI,:));
	goodI = true(size(I)); goodI(nI) = dist(:,1) > -options.geomtol;
	if any(not(goodI))
		warning(sprintf('REFINEMENT CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(goodI))));
		splitedg(splitedg) = goodI;
		newcoords = newcoords(goodI,:);
	end;
end;
ndone = sum(splitedg)/size(edg2fac,1);

I_ = find(splitedg);
I2 = zeros(size(edg2fac,1),1); I2(splitedg) = 1:numel(I_);
deltric = sum(splitedg(tri2edg),2);
deltri = I2(tri2edg);
deltri_ = deltri+size(xy,1);
%deltri contains the splitedg numbers for each element
%deltri_ contains the new node numbers for each element

%UPDATE Nmetric
if size(Nmetric,2) ~= 6
Tmetric = Nmetric(:,7:end); Nmetric = Nmetric(:,1:6); 
Tmetric = [Tmetric; (Tmetric(edg(splitedg,1),:) + Tmetric(edg(splitedg,1),:))/2.];
end;
Nmetric = [Nmetric; metric_avg(edg(splitedg,:),Nmetric,options)]; 
xy = [xy; newcoords];
[xy,Nmetric,keepbndfacs,newfacs,newIDs,delfacc,fac2,fac2q] = mesh_faces(splitedg,I2,edg2fac,fac2edg,fac2edga,edga,bndfac2fac,xy,Nmetric,bndmesh,options); %,tri,fac2tri2
J_ = find(delfacc==2);
delfacc2 = zeros(size(fac,1),1); delfacc2(J_) = 1:numel(J_);

% SPLIT ELEMENTS WITH 1 LONG EDGE
%[tmp,vol] = elem_inv(tri,xy);
I = deltric==1;
if any(I)
    [newtri,newq] = mesh1long(deltri_(I,:),xy,edga,tri2edga(I,:),Nmetric,options);
    oldtri = repmat(find(I),2,1);
    %[tmp1,vol1] = elem_inv(newtri,xy); disp([1 sum(vol1) sum(vol(I))]);
else
    newtri = []; newq = []; oldtri = [];
end;

% SPLIT ELEMENTS WITH 2 LONG EDGES
I = deltric==2;
if any(I)
    [newtri2,q2,oldtri2] = mesh2long(tri(I,:),deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp2,vol2] = elem_inv(newtri2,xy);  disp([2 sum(vol2) sum(vol(I))]);
    newtri = [newtri; newtri2]; newq = [newq; q2]; oldtri = [oldtri; oldtri2]; 
end; 

% SPLIT ELEMENTS WITH 3 LONG EDGES
I = deltric==3;
if any(I)
    [newtri3,q3,xy,Nmetric,oldtri3] = mesh3long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp3,vol3] = elem_inv(newtri3,xy);  disp([3 sum(vol3) sum(vol(I))]);
    newtri = [newtri; newtri3]; newq = [newq; q3]; oldtri = [oldtri; oldtri3];
end; 

%% SPLIT ELEMENTS WITH 4 LONG EDGES
I = deltric==4; 
if any(I)
    [newtri4,q4,xy,Nmetric,oldtri4] = mesh4long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp4,vol4] = elem_inv(newtri4,xy); disp([4 sum(vol4) sum(vol(I))]);
    newtri = [newtri; newtri4]; newq = [newq; q4]; oldtri = [oldtri; oldtri4];
end;

% SPLIT ELEMENTS WITH 5 LONG EDGES
I = deltric==5;
if any(I)
    [newtri5,q5,oldtri5] = mesh5long(deltri_(I,:),xy,edga,tri2edga(I,:),fac2,fac2q,delfacc2,tri2fac(I,:),fac2tri2,find(I),Nmetric,options);
    %[tmp5,vol5] = elem_inv(newtri5,xy); disp([5 sum(vol5) sum(vol(I))]);
    newtri = [newtri; newtri5]; newq = [newq; q5]; oldtri = [oldtri; oldtri5];
end;
% SPLIT ELEMENTS WITH 6 LONG EDGES
I = deltric==6;
if any(I)
    [newtri6,q6] = mesh6long(tri(I,:),deltri_(I,:),xy,Nmetric,options);
  %  [tmp6,vol6] = elem_inv(newtri6,xy); disp([6 sum(vol6) sum(vol(I))]);
    newtri = [newtri; newtri6]; newq = [newq; q6]; 
    oldtri = [oldtri; repmat(find(I),8,1)];
end; 
try; if exist('Tmetric','var'); Nmetric = [Nmetric Tmetric]; end; end %HU


function splitedg_ = splitedg1(splitedg,edgnowL,tri2edg,edg2tri)
edgs = size(splitedg,1);
splitedg_ = false(edgs,1);
splitedgf = zeros(edgs,1); splitedgf(splitedg) = find(splitedg);
edg_ngh = bks_edg2tri2edg(edg2tri,tri2edg); nghI = edg_ngh~=0;
edgnowL(not(splitedg)) = 0;
while any(splitedg)
edg_ngh_ = edg_ngh(splitedg,:);
nI = sum(splitedg);
otherL = zeros(nI,size(edg_ngh,2)); otherL(nghI(splitedg,:)) = edgnowL(edg_ngh_(nghI(splitedg,:)));
otherI = zeros(nI,size(edg_ngh,2)); otherI(nghI(splitedg,:)) = splitedgf(edg_ngh_(nghI(splitedg,:)));
selfL = repmat(edgnowL(splitedg),1,size(edg_ngh,2));
selfI = repmat(splitedgf(splitedg),1,size(edg_ngh,2));
I = false(edgs,1);
I(splitedg) = all(or(selfL > otherL,and(and(selfL~=0, selfL == otherL),selfI>otherI)),2); 
splitedg_(I) = true;
nsplit = edg_ngh(I,:); nsplit = nsplit(nsplit(:)~=0);
splitedg(nsplit(:)) = false; edgnowL(nsplit(:)) = 0;
splitedg(I) = false; splitedgf(I) = 0; 
end;

function [badE,goodE,splitedg_,goodnd] = calc_badE(oldtri,newq,triQ,splitedg,tri2edg,edg2tri)
edg2tri = edg2tri(splitedg,:);
edgQ = inf(size(edg2tri)); edgQ(edg2tri~=0) = triQ(edg2tri(edg2tri~=0));
edgQ_ = zeros(size(splitedg)); edgQ_(splitedg) = min(edgQ,[],2);
newq_ = min(reshape(newq,size(newq,1)/2,2),[],2);
oldtri = oldtri(1:size(oldtri,1)/2);
tri2edg = tri2edg(oldtri,:)'; edg = tri2edg(splitedg(tri2edg));
Ibad = newq_ < edgQ_(edg);
badE = oldtri(Ibad);
if nargout == 1
return;
end;
splitedg_ = splitedg;
splitedg_(edg(Ibad)) = false;
goodnd = splitedg_(splitedg);
goodE = repmat(splitedg_(edg),2,1);




function [abest,bbest] = surfqual(fac,indA,indB,xy,Nmetric,options)
N = size(fac,1); 
qA = elem_qual([fac(:,indA(1,:));fac(:,indA(2,:))],xy,Nmetric,options);
qB = elem_qual([fac(:,indB(1,:));fac(:,indB(2,:))],xy,Nmetric,options);
qAm = min(reshape(qA,N,2),[],2);
qBm = min(reshape(qB,N,2),[],2);
abest = qAm > qBm; bbest= not(abest);



function [xy,Nmetric,keepbndfacs,newfacs,newIDs,delfacc,fac2,fac2q] = mesh_faces(splitedg,I2,edg2fac,fac2edg,fac2edga,edga,bndfac2fac,xy,Nmetric,bndmesh,options) %,tri,fac2tri2
splitfac = false(size(fac2edg,1),1);
splitfac_ = edg2fac(splitedg,:);
splitfac(splitfac_(splitfac_~=0)) = true;
J_ = find(splitfac);
delfacc = sum(splitedg(fac2edg),2); delfacc_ = delfacc(splitfac);
delfac = I2(fac2edg(splitfac,:));
delfac_ = delfac+max(edga(:)); 

N1 = sum(delfacc==1); N2 = sum(delfacc==2); N3 = sum(delfacc==3);
keepbndfacs = true(size(bndfac2fac)); 
bndfac2fac_ = zeros(size(delfacc)); bndfac2fac_(bndfac2fac) = 1:numel(bndfac2fac);
% One split edge per face
if N1~=0

fac2edga1 = fac2edga(delfacc==1,:); %fac1_ = fac(delfacc==1,:);
delfac1_ = delfac_(delfacc_==1,:);
It = delfac(delfacc_==1,:)'~=0; [C,R] = find(It);
fac1 = zeros(N1,4);
fac1(:,1) = delfac1_(R+(C-1)*N1);
fac1(:,[2 3]) = edga(fac2edga1(R+(C-1)*N1),:);
Cp = C+1; Cp(Cp==4) = 1;
fac1(:,4) = edga(fac2edga1(R+(Cp-1)*N1),1);

%bndmesh
It = delfacc(bndfac2fac) == 1; keepbndfacs(It) = false;
bndfac2fac1_ = bndfac2fac_(delfacc==1); bndfac2fac1 = bndfac2fac1_~=0;
newfacs = [fac1(bndfac2fac1,[1 2 4]); fac1(bndfac2fac1,[1 3 4])];
newIDs = repmat(bndmesh.IDs(bndfac2fac1_(bndfac2fac1)),2,1);
else
newfacs = []; newIDs = [];
end;

% Two split edges per face
if N2~=0
delfac2_ = delfac_(delfacc_==2,:);
fac2edga2 = fac2edga(delfacc==2,:);
It = delfac(delfacc_==2,:)'~=0; [C,R] = find(It); [C1,R1] = find(not(It)); C = reshape(C,2,N2);
%cc2 = C1;
%delfac_t2 = delfac_(delfacc_==2,:)';
%fac4 = [edga(fac2edga2(R1+(C(1,:)'-1)*N2),:) edga(fac2edga2(R1+(C(2,:)'-1)*N2),:) C1];
%[all(fac4(or(C1==1,C1==3),1) == fac4(or(C1==1,C1==3),4)) all(fac4(C1==2,2) == fac4(C1==2,3))]
fac2 = zeros(N2,5); fac2(:,[3 4]) = [delfac2_(R1+(C(1,:)'-1)*N2) delfac2_(R1+(C(2,:)'-1)*N2)]; % reshape(delfac_t2(It),2,N2)'; 
It2 = C1==2; I13 = not(It2);
fac2(It2,[3 4]) = fac2(It2,[4 3]);
if any(I13)
fac2(I13,[1 2 5]) = [edga(fac2edga2(R1(I13)+(C(1,I13)'-1)*N2),2) edga(fac2edga2(R1(I13)+(C(2,I13)'-1)*N2),:)]; end;
if any(It2)
fac2(It2,[2 1 5]) = [edga(fac2edga2(R1(It2)+(C(1,It2)'-1)*N2),1) edga(fac2edga2(R1(It2)+(C(2,It2)'-1)*N2),[2 1])]; end;
qualA = elem_qual([fac2(:,[1 2 3]); ...
		 fac2(:,[4 2 3])],xy,Nmetric,options);
qualB = elem_qual([fac2(:,[2 1 4]); ...
		 fac2(:,[3 1 4])],xy,Nmetric,options);
qualAm = min(reshape(qualA,N2,2),[],2);
qualBm = min(reshape(qualB,N2,2),[],2);
fac2q = qualAm > qualBm; fac2nq = not(fac2q);

%[fac2q,fac2nq] = surfqual(fac2,[[1 2 3];[4 2 3]],[[2 1 4]; [3 1 4]],xy,Nmetric,options); 
%tri2 = tri(fac2tri2(delfacc==2),:)';
%It = or(or(tri2 == repmat(fac2(:,1)',4,1),tri2 == repmat(fac2(:,2)',4,1)),tri2 == repmat(fac2(:,5)',4,1));
%loneND = tri2(not(It));
%[tmp,vol] = elem_inv([fac2(:,[1 2 5]) loneND],xy); [tmp cc2]

%bndmesh
It = delfacc(bndfac2fac) == 2; keepbndfacs(It) = false;
bndfac2fac2_ = bndfac2fac_(delfacc==2); bndfac2fac2 = bndfac2fac2_~=0;
ItA_ = and(bndfac2fac2,fac2q ); ItA = ItA_(bndfac2fac2);
ItB_ = and(bndfac2fac2,fac2nq); ItB = ItB_(bndfac2fac2);
newfacs2 = zeros(sum(bndfac2fac2),9);
newfacs2(ItA,[1 4 7]) = fac2(ItA_,[1 2 3]);
newfacs2(ItA,[2 5 8]) = fac2(ItA_,[4 2 3]);
newfacs2(ItB,[1 4 7]) = fac2(ItB_,[2 1 4]);
newfacs2(ItB,[2 5 8]) = fac2(ItB_,[3 1 4]);
newfacs2(:,[3 6 9]) = fac2(bndfac2fac2,[3 4 5]);
newfacs = [newfacs; reshape(newfacs2,sum(It)*3,3)];
newIDs = [newIDs; repmat(bndmesh.IDs(bndfac2fac2_(bndfac2fac2)),3,1)];
else
fac2 = []; fac2q = []; cc2 = [];
end;
if N3~=0
delfac3_ = delfac_(delfacc_==3,:);
fac2edga3 = fac2edga(delfacc==3,:);
fac3 = zeros(N3,6);
fac3(:,[1 2 3]) = delfac3_;
fac3(:,[4 5]) = edga(fac2edga3(:,3),[2 1]);
fac3(:,6) = edga(fac2edga3(:,1),1);
%21=32 opposite of edga1
%31=12 opposite of edga2
%11=22 opposite of edga3

%bndmesh
It = delfacc(bndfac2fac) == 3; keepbndfacs(It) = false;
bndfac2fac3_ = bndfac2fac_(delfacc==3); bndfac2fac3 = bndfac2fac3_~=0;
newfacs = [newfacs; fac3(bndfac2fac3,[1 2 3]); fac3(bndfac2fac3,[2 3 4]); fac3(bndfac2fac3,[3 1 5]); fac3(bndfac2fac3,[1 2 6])];
newIDs = [newIDs; repmat(bndmesh.IDs(bndfac2fac3_(bndfac2fac3)),4,1)];
end;
%bndmesh.fac = [bndmesh.fac(keepbndfacs,:); sort(newfacs,2)];
%bndmesh.IDs = [bndmesh.IDs(keepbndfacs); newIDs];



function [newtri1,q1] = mesh1long(deltri1_,xy,edga,tri2edga1,Nmetric,options)
%deltri1_ = deltri_(I,:);
%tri2edga1 = tri2edga(I,:); 
tzr = max(edga(:)); N1= size(tri2edga1,1);
Ca1  = [1 3 5 6 9 11]';
deltri1_n = deltri1_';
I1 = deltri1_n~=tzr;
[C,R] = find(I1);
spltnd = deltri1_n(I1);
edgaC1 = tri2edga1(R+(Ca1(C)-1)*N1);
thmsk = [2 1 4 3 6 5]';
edgaC2 = tri2edga1(R+(Ca1(thmsk(C))-1)*N1); 
nodes = [edga(edgaC1,:) edga(edgaC2,:)];
flp = or(C==1,C==2); nodes(flp,[3 4]) = nodes(flp,[4 3]);
newtri1 = [[spltnd nodes(:,1) nodes(:,[3 4])]; ...
	  [spltnd nodes(:,2) nodes(:,[4 3])]];
q1 = elem_qual(newtri1,xy,Nmetric,options,1);

%[tmp1,vol1] = elem_inv(newtri1,xy);
%[tmp,vol] = elem_inv(tri(I,:),xy);
%[sum(vol) sum(vol1)]



function [newtri2,q2,oldtri] = mesh2long(tri2, deltri2_,xy,edga,tri2edga2,fac2,fac2q,delfacc2,tri2fac2,fac2tri2,fI2,Nmetric,options)
%deltri2_ = deltri_(I,:); 
%tri2edga2 = tri2edga(I,:); 
%tri2=tri(I,:); 
tzr = max(edga(:)); N2 = size(tri2,1);
Ca1  = [1 3 5 6 9 11]';
I_ = [all(deltri2_(:,[1 2]) ~= tzr,2) all(deltri2_(:,[3 4]) ~= tzr,2) all(deltri2_(:,[5 6]) ~= tzr,2)];
EA = any(I_,2); nEA = not(EA);

%CASE 1: two opposing split edges
if any(EA)
deltri2n_ = deltri2_(EA,:); N2 = size(deltri2n_,1);
tri2edga2n = tri2edga2(EA,:);
[C,R] = find(I_(EA,:)');
C2 = 2*C;
spltnds = [deltri2n_(R+(C2-1)*N2) deltri2n_(R+(C2-2)*N2)];
crnrs = [edga(tri2edga2n(R+(Ca1(C2)-1)*N2),:) edga(tri2edga2n(R+(Ca1(C2-1)-1)*N2),:)];
flp = C==1; spltnds(flp,:) = spltnds(flp,[2 1]);
newtri2 = [[spltnds crnrs(:,1) crnrs(:,3)]; ...
	  [spltnds crnrs(:,4) crnrs(:,1)]; ...
  	  [spltnds crnrs(:,3) crnrs(:,2)]; ...
	  [spltnds crnrs(:,2) crnrs(:,4)]];
%[tmp2,vol2] = elem_inv(newtri2,xy); 
%[tmp,vol] = elem_inv(tri(I,:),xy); vol=vol(EA);
%[sum(vol) sum(vol2)]
q2 = elem_qual(newtri2,xy,Nmetric,options);
oldtri = repmat(fI2(EA),4,1);
else
newtri2 = []; q2 = []; oldtri = [];
end;%ifany EA

%CASE 2: two opposing edges not split
if any(nEA) 
deltri2n_ = deltri2_(nEA,:); 
tri2edga2n = tri2edga2(nEA,:);
tri2n = tri2(nEA,:); N2 = size(tri2n,1);

tri2fac2n = tri2fac2(nEA,:); fact = reshape(delfacc2(tri2fac2n)',4,N2); %reshape relevant for N2==1
[Cn,R] = find(fact); 
trinEA = fI2(nEA); fact_ = tri2fac2n'; fact_ = fact_(fact~=0); trinEA_ = fac2tri2(fact_);
fact = fact(fact~=0); qAm = fac2q(fact); fac = fac2(fact,1:4); shrnd = fac2(fact,5); % CC2 = cc2(fact);
flp_ = trinEA_ == trinEA; 
%fact2 = fact(:,[2 1 4 3]);
loneCR = zeros(N2,1); 
I_ = Cn==3; loneCR(I_) = tri2n(I_,1);
I_ = Cn==4; loneCR(I_) = tri2n(I_,2);
I_ = Cn==1; loneCR(I_) = tri2n(I_,3);
I_ = Cn==2; loneCR(I_) = tri2n(I_,4);

deltri2n_t = deltri2n_'; I_ = deltri2n_t~=tzr; 
spltnds = reshape(deltri2n_t(I_),2,N2)';
[C,R] = find(I_); R=R(1:2:end); C1=C(1:2:end); C2=C(2:2:end); Call=C1+(C2-1)*4;

%CallN = Call+(CC2-1)*24;
%thmsk1=false(72,1); thmsk1([41 45 46]) = true;
%thmsk2=false(72,1); thmsk2([10 18 19 20 24 62 65 66 69 71]) = true;
%flp = or(and(thmsk1(CallN),flp_),and(thmsk2(CallN),not(flp_))); fac(flp,:) = fac(flp,[2 1 4 3]);
thmsk=false(24,1); thmsk([17 18 21 22]) = true;
flp = thmsk(Call); spltnds(flp,:) = spltnds(flp,[2 1]);
flp = fac(:,4) == spltnds(:,1);
fac(flp,:) = fac(flp,[2 1 4 3]);
qAm(flp) = not(qAm(flp)); nqAm = not(qAm);


%%[all(fact==fac,2) all(fact2==fac,2)]
%%CallN_ = CallN(flp_);
%%unique(CallN_(all(fact(flp_,:)==fac(flp_,:),2)))
%%unique(CallN_(not(all(fact(flp_,:)==fac(flp_,:),2))))
%
%crnrs = [edga(tri2edga2n(R+(Ca1(C1)-1)*N2),:) ...
	%edga(tri2edga2n(R+(Ca1(C2)-1)*N2),:)];
%%put shrnd in crnrs(:,[2 3])
%thmsk=false(24,1); thmsk([9 10 18 19 21 24]) = true;
%flp = thmsk(Call); crnrs(flp,[1 2]) = crnrs(flp,[2 1]);
%thmsk=false(24,1); thmsk([9 13 18 19 21 23]) = true;
%flp = thmsk(Call); crnrs(flp,[3 4]) = crnrs(flp,[4 3]);
%shrnd = crnrs(:,2);
%loneCR = zeros(N2,1);
%I_ = any([Call==10 Call==22 Call==23],2); loneCR(I_) = tri2n(I_,1);
%I_ = any([Call==14 Call==18 Call==20],2); loneCR(I_) = tri2n(I_,2);
%I_ = any([Call==9  Call==17 Call==19],2); loneCR(I_) = tri2n(I_,3);
%I_ = any([Call==13 Call==21 Call==24],2); loneCR(I_) = tri2n(I_,4);
%%crnrs_ = sort(crnrs'); shrnd_=crnrs_([false(1,N2);  crnrs_(1:end-1,:)==crnrs_(2:end,:)]);
%%I = [any(repmat(tri2n(:,1),1,3) == crnrs(:,[1 2 4]),2) ...
%%     any(repmat(tri2n(:,2),1,3) == crnrs(:,[1 2 4]),2) ...
%%     any(repmat(tri2n(:,3),1,3) == crnrs(:,[1 2 4]),2) ...
%%     any(repmat(tri2n(:,4),1,3) == crnrs(:,[1 2 4]),2)]';
%%tri2n = tri2n'; loneCR_ = tri2n(not(I));
%fac = [crnrs(:,[1 4]) spltnds];
%thmsk=false(24,1); thmsk([17 18 21 22]) = true;
%flp = thmsk(Call); fac(flp,:) = fac(flp,[2 1 4 3]);
%%[tmpa,areaa] = elem_inv([[fac(:,[1 2 3]); fac(:,[4 2 3])]],xy); %areaa = reshape(areaa,N2,2);
%%[tmpb,areab] = elem_inv([[fac(:,[2 1 4]); fac(:,[3 1 4])]],xy); %areab = reshape(areab,N2,2);
%%[tmp,area] = elem_inv([crnrs(:,[1 4]) shrnd],xy); 
%%[area./sum(areaa,2)*3/4 area./sum(areab,2)*3/4]
		 %
%%qualA = elem_qual([fac(:,[1 2 3]); ...
		 %%fac(:,[4 2 3])],xy,Nmetric,options);
%%qualB = elem_qual([fac(:,[2 1 4]); ...
		 %%fac(:,[3 1 4])],xy,Nmetric,options);
%%qualAm = min(reshape(qualA,N2,2),[],2);
%%qualBm = min(reshape(qualB,N2,2),[],2);
%%qAm = qualAm > qualBm; nqAm = not(qAm);
%[qAm,nqAm] = surfqual(fac,[[1 2 3];[4 2 3]],[[2 1 4]; [3 1 4]],xy,Nmetric,options); 
for i=1:2
 if i==1
 pos=qAm; j2=2; j3=3;
 else
 pos=nqAm; j2=1; j3=4;
 end;
 if not(any(pos))
 continue;
 end;
 newtri2_ = [[loneCR(pos) fac(pos,[4 3]) shrnd(pos)]; ...
 	     [loneCR(pos) fac(pos,[j2 3 4])]; ...
	     [loneCR(pos) fac(pos,[2 1 j3])]];
 %[tmp1,vol1] = elem_inv(newtri2_,xy); tmp1 = reshape(tmp1,sum(pos),3); vol1 = reshape(vol1,sum(pos),3); 
 %[tmp,vol] = elem_inv(tri2n,xy); %Call_ = Call(pos); 
 %[sum(vol1(:)) sum(vol(pos))]
 newtri2 = [newtri2; newtri2_]; 
 q2 = [q2; elem_qual(newtri2_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(trinEA(pos),3,1)];
end;%fori
end;%ifany nEA



function [newtri3,q3,xy,Nmetric,oldtri] = mesh3long(deltri3_,xy,edga,tri2edga3,fac2,fac2q,delfacc2,tri2fac3,fac2tri2,fI3,Nmetric,options)
%deltri3_ = deltri_(I,:);
%tri2edga3 = tri2edga(I,:);
tzr = max(edga(:));
edgs = size(edga,1);
Ca1  = [1 3 5 6 9 11]'; %find([true false true false true true false false true false true false]) (gen_books3D)

EA = all([sum(deltri3_(:,[1 2])~=tzr,2) ...
          sum(deltri3_(:,[3 4])~=tzr,2) ...
          sum(deltri3_(:,[5 6])~=tzr,2)]==1,2); nEA = not(EA);
          
% CASE 1+2
if any(EA) %no opposing split edges
deltri3n_ = deltri3_(EA,:);
tri2edga3n = tri2edga3(EA,:); 
fI3n = fI3(EA); N3 = size(deltri3n_,1);
[C1,R1] = find(deltri3_(EA,[1 2])' ~= tzr);
[C2,R2] = find(deltri3_(EA,[3 4])' ~= tzr);
[C3,R3] = find(deltri3_(EA,[5 6])' ~= tzr); 
Call = C1+(C2-1)*2+(C3-1)*4; 
C2=C2+2; C3=C3+4; C1p = C1+1; C1p(C1p==3) = 1;
thmsk1 = false(8,1); thmsk1([1 4 6 7]) = true;
thmsk2 = not(thmsk1); I1 = thmsk1(Call); I2=thmsk2(Call);
if any(I1)
R1 = R1(I1); N1 = numel(R1);
thmsk = [2 0 0 1 0 2 1 0]';
CE1 = thmsk(Call(I1));
loneCR1 = edga(tri2edga3n(R1+(Ca1(C1p(I1))-1)*N3)+(CE1-1)*edgs);
%%brute find
%C2p = C2+1; C2p(C2p==5) = 3;
%C3p = C3+1; C3p(C3p==7) = 5;
%nodes_ = sort([edga(tri2edga3n(R1+(Ca1(C1p(I1))-1)*N3),:) ...
	      %edga(tri2edga3n(R1+(Ca1(C2p(I1))-1)*N3),:) ...
	      %edga(tri2edga3n(R1+(Ca1(C3p(I1))-1)*N3),:)],2)';
%loneCR_ = nodes_([nodes_(2:end,:)==nodes_(1:end-1,:); false(1,N1)]); loneCR_ = loneCR_(repmat([false; true],N1,1));
%[loneCR_ loneCR1 Call(I1)]

ncrnrs = [deltri3n_(R1+(C1(I1)-1)*N3) deltri3n_(R1+(C2(I1)-1)*N3) deltri3n_(R1+(C3(I1)-1)*N3)];
thmska = [1 0 0 2 0 1 2 0]';
thmskb = [2 0 0 1 0 2 1 0]';
thmskc = [1 0 0 2 0 2 1 0]';
CE1a = thmska(Call(I1));
CE1b = thmskb(Call(I1));
CE1c = thmskc(Call(I1));
crnrs = [edga(tri2edga3n(R1+(Ca1(C2(I1))-1)*N3)+(CE1c-1)*edgs) ...
	edga(tri2edga3n(R1+(Ca1(C1(I1))-1)*N3)+(CE1b-1)*edgs) ...
	edga(tri2edga3n(R1+(Ca1(C1(I1))-1)*N3)+(CE1a-1)*edgs)];
newtri3 = [[loneCR1 ncrnrs]; ...
	  [loneCR1 ncrnrs(:,[2 1]) crnrs(:,3)]; ...
	  [loneCR1 ncrnrs(:,[3 2]) crnrs(:,1)]; ...
	  [loneCR1 ncrnrs(:,[1 3]) crnrs(:,2)]];
q3 = elem_qual(newtri3,xy,Nmetric,options);
oldtri = repmat(fI3n(I1),4,1);
%[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA);
%[tmp_,vol_] = elem_inv(newtri3,xy);
%repmat(vol(I1),4,1)./vol_
%[sum(vol(I1)) sum(vol_)]
else
newtri3 = []; q3 = []; oldtri = [];
end; %ifany I1

if any(I2)
R2 = R2(I2); N2 = nnz(I2);
ItF = [true(N2,1); false(N2,1); false(N2,1)]; 
ItM = [false(N2,1); true(N2,1); false(N2,1)]; 
ItL = [false(N2,1); false(N2,1); true(N2,1)];

%ncrnrs = [deltri3n_(R2+(C1(I2)-1)*N3) deltri3n_(R2+(C2(I2)-1)*N3) deltri3n_(R2+(C3(I2)-1)*N3)];
%thmsk = [0 1 2 0 1 0 0 2]';
%CE2 = thmsk(Call(I2));
%loneCR2 = edga(tri2edga3n(R2+(Ca1(C1(I2))-1)*N3)+(CE2-1)*edgs);
%%nodes_ = sort([edga(tri2edga3n(R2+(Ca1(C1(I2))-1)*N3),:) ...
%%	      edga(tri2edga3n(R2+(Ca1(C2(I2))-1)*N3),:) ...
%%	      edga(tri2edga3n(R2+(Ca1(C3(I2))-1)*N3),:)],2)';
%%loneCR_ = nodes_([nodes_(2:end,:)==nodes_(1:end-1,:); false(1,N2)]); loneCR_ = loneCR_(repmat([false; true],N2,1));
%%[loneCR_ loneCR2 Call(I2)]

tri2fac3t = tri2fac3(EA,:)'; tri2fac3t = tri2fac3t(:,I2); fI3nn = fI3n(I2); fact = delfacc2(tri2fac3t);
fact_ = tri2fac3t; fact_ = reshape(tri2fac3t(fact~=0),3,N2)';
trinEA_ = reshape(fac2tri2(fact_),N2,3); %reshape, when N2==1
fact = reshape(fact(fact~=0),3,N2)'; qAm = fac2q(fact(:)); fac = fac2(fact(:),1:4); loneCR2 = fac2(fact(:,1),5);
flp_ = trinEA_ == repmat(fI3nn,1,3);
fac(flp_(:),:) = fac(flp_(:),[2 1 4 3]); qAm(flp_(:)) = not(qAm(flp_(:)));
thmsk = false(8,1); thmsk([2 5]) = true;
It_ = thmsk(Call(I2)); ItF_ = [It_; false(N2*2,1)]; ItM_ = [false(N2,1); It_; false(N2,1)];
tmp = fac(ItF_,:); fac(ItF_,:) = fac(ItM_,:); fac(ItM_,:) = tmp;
tmp = qAm(ItF_);   qAm(ItF_)   = qAm(ItM_);   qAm(ItM_)   = tmp;
ncrnrs = [fac(ItF,4) fac(ItM,4) fac(ItL,4)];
crnrs  = [fac(ItF,2) fac(ItM,2) fac(ItL,2)];
nqAm = qAm; qAm = not(nqAm); fac = fac(:,[4 3 2 1]); %could be avoided


%thmska = [0 2 1 0 2 0 0 1]';
%thmskb = [0 1 2 0 1 0 0 2]';
%CE1a = thmska(Call(I2));
%CE1b = thmskb(Call(I2));
%crnrs = [edga(tri2edga3n(R2+(Ca1(C1(I2))-1)*N3)+(CE1a-1)*edgs) ...
	%edga(tri2edga3n(R2+(Ca1(C1p(I2))-1)*N3)+(CE1b-1)*edgs) ...
	%edga(tri2edga3n(R2+(Ca1(C1p(I2))-1)*N3)+(CE1a-1)*edgs)];
%%[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA);
%%[tmp_,vol_] = elem_inv([loneCR2 ncrnrs],xy);
%%[sum(vol(I2)) sum(vol_)*8]
%%III=4; figure(); hold on; 
%%plot3(xy(ncrnrs(III,[1 2 3 1]),1),xy(ncrnrs(III,[1 2 3 1]),2),xy(ncrnrs(III,[1 2 3 1]),3),'-b.'); 
%%plot3(xy(crnrs(III,[1 2 3 1]),1),xy(crnrs(III,[1 2 3 1]),2),xy(crnrs(III,[1 2 3 1]),3),'-k.');
%%plot3([xy(ncrnrs(III,2),1) xy(crnrs(III,2),1)],[xy(ncrnrs(III,2),2) xy(crnrs(III,2),2)],[xy(ncrnrs(III,2),3) xy(crnrs(III,2),3)],'g-');
%%plot3([xy(ncrnrs(III,1),1) xy(crnrs(III,1),1)],[xy(ncrnrs(III,1),2) xy(crnrs(III,1),2)],[xy(ncrnrs(III,1),3) xy(crnrs(III,1),3)],'r-'); hold off;
%
newtri3 = [newtri3; [loneCR2 ncrnrs]];
q3 = [q3; elem_qual([loneCR2 ncrnrs],xy,Nmetric,options)];
oldtri = [oldtri; fI3nn];
%
%fac1 = [ncrnrs(:,[1 2]) crnrs(:,[1 2])];
%fac2 = [ncrnrs(:,[2 3]) crnrs(:,[2 3])];
%fac3 = [ncrnrs(:,[3 1]) crnrs(:,[3 1])];
%fac = [fac1;fac2;fac3];
%%[tmp1_,vol1_] = elem_inv([fac1(:,[1 2 3]);fac1(:,[1 3 4])],xy); vol1_ = sum(reshape(vol1_,N2,2),2);
%%[tmp2_,vol2_] = elem_inv([fac2(:,[1 2 3]);fac2(:,[1 3 4])],xy); vol2_ = sum(reshape(vol2_,N2,2),2);
%%[tmp3_,vol3_] = elem_inv([fac3(:,[1 2 3]);fac3(:,[1 3 4])],xy); vol3_ = sum(reshape(vol3_,N2,2),2);
%%[tmpa,vol1a] = elem_inv([crnrs(:,[1 2]) loneCR2],xy);
%%[tmpb,vol1b] = elem_inv([ncrnrs(:,[1 2]) loneCR2],xy);
%%[tmpa,vol2a] = elem_inv([crnrs(:,[2 3]) loneCR2],xy);
%%[tmpb,vol2b] = elem_inv([ncrnrs(:,[2 3]) loneCR2],xy);
%%[tmpa,vol3a] = elem_inv([crnrs(:,[3 1]) loneCR2],xy);
%%[tmpb,vol3b] = elem_inv([ncrnrs(:,[3 1]) loneCR2],xy);
%%[sum(vol1a) sum(vol1b)*4 sum(vol1_)*4/3]
%%[sum(vol2a) sum(vol2b)*4 sum(vol2_)*4/3]
%%[sum(vol3a) sum(vol3b)*4 sum(vol3_)*4/3]
%
%%qualA = elem_qual([fac(:,[2 1 4]); ...
		 %%fac(:,[3 1 4])],xy,Nmetric,options);
%%qualB = elem_qual([fac(:,[1 2 3]); ...
		 %%fac(:,[4 2 3])],xy,Nmetric,options);
%%qualAm = min(reshape(qualA,N2*3,2),[],2);
%%qualBm = min(reshape(qualB,N2*3,2),[],2);
%%qAm = qualAm > qualBm; nqAm = not(qAm);
%[qAm,nqAm] = surfqual(fac,[[2 1 4];[3 1 4]],[[1 2 3]; [4 2 3]],xy,Nmetric,options); 

%no symmetry (75 % of cases)
It1 = all([qAm(ItF) ,qAm(ItM), nqAm(ItL)],2); %1 to 4
It2 = all([nqAm(ItF),nqAm(ItM),qAm(ItL) ],2); %2 to 3
It3 = all([qAm(ItM) ,qAm(ItL), nqAm(ItF)],2); %1 to 4
It4 = all([nqAm(ItM),nqAm(ItL),qAm(ItF) ],2); %2 to 3
It5 = all([qAm(ItL) ,qAm(ItF), nqAm(ItM)],2); %1 to 4
It6 = all([nqAm(ItL),nqAm(ItF),qAm(ItM) ],2); %2 to 3
for i=1:6
switch i
case 1
pos_=It1; ItFMLa = ItM; ItFMLb = ItF; k1=1; k2=2; j1=1; j4=4;
case 2
pos_=It2; ItFMLa = ItF; ItFMLb = ItM; k1=2; k2=1; j1=3; j4=2;
case 3
pos_=It3; ItFMLa = ItL; ItFMLb = ItM; k1=1; k2=2; j1=1; j4=4;
case 4
pos_=It4; ItFMLa = ItM; ItFMLb = ItL; k1=2; k2=1; j1=3; j4=2;
case 5
pos_=It5; ItFMLa = ItF; ItFMLb = ItL; k1=1; k2=2; j1=1; j4=4;
case 6
pos_=It6; ItFMLa = ItL; ItFMLb = ItF; k1=2; k2=1; j1=3; j4=2;
end;
if not(any(pos_))
	continue;
end;
pos = repmat(pos_,3,1); k3=k1+2; k4=k2+2;
ItFMLa_ = and(pos,ItFMLa); ItFMLb_ = and(pos,ItFMLb);
newtri3_ = [[fac(ItFMLa_,[j1 j4 k2]) fac(ItFMLb_,k1)]; ...
	   [fac(ItFMLb_,[j4 j1 k3]) fac(ItFMLa_,k4)]; ...
	   [fac(ItFMLb_,[j1 j4]) fac(ItFMLa_,[k1 k4])]];
%[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA); vol = vol(I2);
%[tmp_,vol_] = elem_inv(newtri3_,xy); 
%[i any(tmp_) sum(vol(pos(ItF))) sum(vol_)*8/7]
newtri3 = [newtri3; newtri3_];
q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
oldtri = [oldtri; repmat(fI3nn(pos_),3,1)];
end; %fori

%rotational symmetry (impossible case)
ItA = all(reshape(qAm,N2,3),2); %1 to 4
ItB = all(reshape(nqAm,N2,3),2); %2 to 3
if any(or(ItA,ItB))
for i=1:2
if i==1
pos_ = ItA; i1 = 1; i2=4;
else
pos_ = ItB; i1 = 2; i2=3;
end;
if not(any(pos_))
	continue;
end;
In = [ncrnrs(pos_,:) crnrs(pos_,:)];
nxy = (xy(In(:,1),:)+xy(In(:,2),:)+xy(In(:,3),:) ...
    +  xy(In(:,4),:)+xy(In(:,5),:)+xy(In(:,6),:))/6.;
nNmetric1 =  metric_avg(In(:,[1 2 3]),Nmetric,options);
nNmetric2 =  metric_avg(In(:,[4 5 6]),Nmetric,options);
edg_ = (1:size(nxy,1))';
nNmetric = metric_avg([edg_ edg_+(edg_(end))],[nNmetric1; nNmetric2],options);
pos = repmat(pos_,3,1); %crnrs(pos,:) ncrnrs(pos,:)
ItF_ = and(pos,ItF); ItM_ = and(pos,ItM); ItL_ = and(pos,ItL);
edg_ = edg_ + size(xy,1);
newtri3_ = [[edg_ fac(ItF_,1) fac(ItL_,1) fac(ItM_,1)]; ... 
	   [edg_ fac(ItF_,3) fac(ItM_,3) fac(ItL_,3)]; ... 
	   [edg_ fac(ItF_,[1 2 i2])]; ...
	   [edg_ fac(ItF_,[i1 4 3])]; ...
	   [edg_ fac(ItM_,[1 2 i2])]; ...
	   [edg_ fac(ItM_,[i1 4 3])]; ...
	   [edg_ fac(ItL_,[1 2 i2])]; ...
	   [edg_ fac(ItL_,[i1 4 3])]];
xy_ = [xy; nxy];
%[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA); vol = vol(I2);
%[tmp_,vol_] = elem_inv(newtri3_,xy_); vol_ = %sum(reshape(vol_,sum(ItF_),8),2);
%[i any(tmp_) sum(vol(pos(ItF))) sum(vol_)*8/7]
xy = [xy; nxy];
Nmetric = [Nmetric; nNmetric];
newtri3 = [newtri3; newtri3_];
q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
oldtri = [oldtri; repmat(fI3nn(pos_),8,1)];
end; %fori
end; %ifany
end; %ifany I2
else
newtri3 = []; q3 =[]; oldtri = [];
end; %ifany EA

%[tmp,vol] =elem_inv(tri(I,:),xy); vol=vol(EA);
%[tmp1,vol1] =elem_inv(newtri3,xy);
%[sum(vol) sum(vol1)]

% CASE 3
if any(nEA) %two opposing split edges 
deltri3n_ = deltri3_(nEA,:);
tri2edga3n = tri2edga3(nEA,:); 
fI3n = fI3(nEA); N3 = size(deltri3n_,1);
[C,R] = find([all(deltri3n_(:,[1 2])~=tzr,2) all(deltri3n_(:,[3 4])~=tzr,2) all(deltri3n_(:,[5 6])~=tzr,2)]');
C = 2*C;
I1 = R+(C-1)*N3; I2 = R+(C-2)*N3;
oppnds = [deltri3n_(I1) deltri3n_(I2)];
deltri3n_([I1;I2])=tzr; deltri3n_ = deltri3n_';
I_ = deltri3n_~=tzr;
[Cs,Rs] = find(I_);
Call = Cs+((C/2)-1)*6;
ItF = [true(N3,1); false(N3,1)]; ItL = not(ItF);

tri2fac3t = tri2fac3(nEA,:)'; fact = delfacc2(tri2fac3t);
fact_ = tri2fac3t; fact_ = reshape(tri2fac3t(fact~=0),2,N3)';
trinEA_ = reshape(fac2tri2(fact_),N3,2); %reshape, when N2==1
fact = reshape(fact(fact~=0),2,N3)'; qsm1 = fac2q(fact(:)); fac = fac2(fact(:),1:4);
flp_ = trinEA_ == repmat(fI3n,1,2);
fac(flp_(:),:) = fac(flp_(:),[2 1 4 3]); qsm1(flp_(:)) = not(qsm1(flp_(:)));

%destroy information by %enforcing fac(ItF,3) == fac(ItL,3)
thmsk = false(18,1); thmsk([5 6 7 8 15 16]) = true;
flp = repmat(thmsk(Call),2,1); fac(flp,:) = fac(flp,[2 1 4 3]);
qsm1(flp) = not(qsm1(flp));
nqsm1 = not(qsm1); fac = fac(:,[3 4 1 2]); %nqsm1 = qsm1; qsm1 = not(nqsm1);%could be avoided


%%thmsk = [2 1 4 3 6 5]';
%%crnrs = [edga(tri2edga3n(Rs+(Ca1(Cs)-1)*N3),:) edga(tri2edga3n(Rs+(Ca1(thmsk(Cs))-1)*N3),:)]; %34 connects the cleaner faces
%shrnd = deltri3n_(I_);
%deltri3n_ = deltri3_(nEA,:);
%I_ = [all(deltri3n_(:,[1 2])==tzr,2) all(deltri3n_(:,[3 4])==tzr,2) all(deltri3n_(:,[5 6])==tzr,2)];
%[C_,R_] = find(I_'); C_ = 2*C_;
%crnrs2 = [edga(tri2edga3n(R_+(Ca1(C_)-1)*N3),:) edga(tri2edga3n(R_+(Ca1(C_-1)-1)*N3),:)]; %12 and 34 are opposing edges not refined
%thmsk = false(18,1); thmsk([4 6 8 11 14 15]) = true; 
%flp = thmsk(Call); crnrs2(flp,[1 2]) = crnrs2(flp,[2 1]);
%%any(any(crnrs2(:,[1 1]) == crnrs(:,[1 2]),2))
%thmsk = false(18,1); thmsk([4 5 7 12 14 15]) = true; 
%flp = thmsk(Call); crnrs2(flp,[3 4]) = crnrs2(flp,[4 3]);
%%any(any(crnrs2(:,[3 3]) == crnrs(:,[1 2]),2))
%thmsk = false(18,1); thmsk([4 6 8 12 14 16]) = true; 
%flp = thmsk(Call); crnrs2(flp,:) = crnrs2(flp,[3 4 1 2]);
%
%fac1 = [shrnd oppnds(:,1) crnrs2(:,[2 1])];
%fac2 = [shrnd oppnds(:,2) crnrs2(:,[4 3])];
%%[tmp1,vol1_] = elem_inv([fac1(:,[1 2 3]); fac1(:,[2 3 4])],xy); vol1_ = sum(reshape(vol1_,N3,2),2);
%%[tmp2,vol2_] = elem_inv([fac2(:,[1 2 3]); fac2(:,[2 3 4])],xy); vol2_ = sum(reshape(vol2_,N3,2),2);
%%[tmp1,vol1] = elem_inv([crnrs(:,[1 2]) fac1(:,4)],xy);
%%[tmp2,vol2] = elem_inv([crnrs(:,[1 2]) fac2(:,4)],xy);
%%[abs(vol1./vol1_*3/4-1)<1e-12 abs(vol2./vol2_*3/4-1)<1e-12 Call]
%fac = [fac1; fac2];
%
%%quals1 = elem_qual([fac(:,[2 1 4]); ...
		  %%fac(:,[3 1 4])],xy,Nmetric,options);
%%quals2 = elem_qual([fac(:,[1 2 3]); ...
		  %%fac(:,[4 2 3])],xy,Nmetric,options);
%%quals1m = min(reshape(quals1,N3*2,2),[],2);
%%quals2m = min(reshape(quals2,N3*2,2),[],2);
%%qsm1 = quals1m > quals2m; nqsm1 = not(qsm1);
%[qsm1,nqsm1] = surfqual(fac,[[2 1 4];[3 1 4]],[[1 2 3]; [4 2 3]],xy,Nmetric,options); 

%first 50% (rotational symmetry)
It1 = all(reshape(qsm1,N3,2),2); %all to shrnd
It2 = all(reshape(nqsm1,N3,2),2); %none to shrnd 
thmsk = false(18,1); thmsk([3 4 11 12 13 14]) = true; flp = thmsk(Call); nflp = not(flp);
for i=1:4
switch i
case 1
 pos_ = and(flp,It1); m1=1; m4=4; j2=2; j3=3; k2=2; k4=4;
case 2
 pos_ = and(nflp,It1); m1=4; m4=1; j2=2; j3=3; k2=4; k4=2;
case 3
 pos_ = and(flp,It2); m1=2; m4=3; j2=4; j3=1; k4=4; k2=2;
case 4
 pos_ = and(nflp,It2);  m1=3; m4=2; j2=4; j3=1; k4=2; k2=4;
end;
if not(any(pos_))
continue;
end;
 pos = repmat(pos_,2,1);
 ItF_ = and(pos,ItF); ItL_ = and(pos,ItL);
 TriT1 = [[fac(ItF_,2) fac(ItL_,[m4 m1 j3])]; ...
         [fac(ItL_,2) fac(ItF_,[m4 m1 j3])]; ...
         [fac(ItF_,2) fac(ItL_,[m1 m4 j2])]; ...
         [fac(ItL_,2) fac(ItF_,[m1 m4 j2])]; ...
         [fac(ItL_,[2 4]) fac(ItF_,[k2 k4])]];
 oldtri1 = repmat(fI3n(pos_,:),5,1);
 if i==1 || i==2
 N3_ = sum(ItF_);
 TriT2 = [[fac(ItF_,2) fac(ItL_,[m1 3 m4])]; ...
         [fac(ItL_,2) fac(ItF_,[m1 3 m4])]; ...
         [fac(ItF_,4) fac(ItL_,[m4 2 m1])]; ...
         [fac(ItL_,4) fac(ItF_,[m4 2 m1])]];
 oldtri2 = repmat(fI3n(pos_,:),4,1);
 qual1 = elem_qual(TriT1,xy,Nmetric,options);
 qual2 = elem_qual(TriT2,xy,Nmetric,options);
 qual1m = min(reshape(qual1,N3_,5),[],2);
 qual2m = min(reshape(qual2,N3_,4),[],2);
 option1 = qual1m > qual2m; option2 = not(option1);
 option1 = repmat(option1,5,1);
 option2 = repmat(option2,4,1);
 newtri3_ = [TriT1(option1,:); TriT2(option2,:)];
 oldtri_ = [oldtri1(option1); oldtri2(option2)];
 q3 = [q3; qual1(option1,:); qual2(option2,:)];
 else
 	newtri3_ = TriT1;
 	oldtri_ = oldtri1;
 	q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
 end;
 %[tmp1,vol1] =elem_inv(newtri3_,xy); vol1=reshape(vol1,sum(ItF_),4); tmp1=reshape(tmp1,sum(ItF_),4);
 %[tmp,vol] =elem_inv(tri(I,:),xy); vol=vol(nEA); %plot(tmp1)
 %[sum(vol(pos(ItF))) sum(vol1(:))]
  %[vol(pos(ItF)) sum(vol1,2)]
 % [Call(pos(ItF)) any(tmp1,2)]
 %Call_ = Call(pos(ItF)); unique(Call_(not(any(tmp1,2))))
 newtri3 = [newtri3; newtri3_];
 oldtri = [oldtri; oldtri_];
end; %fori

%no symmetry
ItA = and(qsm1(ItF),nqsm1(ItL));
ItB = and(nqsm1(ItF),qsm1(ItL));
if any([ItA(:);ItB(:)])
for i=1:4
switch i
case 1
 pos_ = and(flp,ItA); ItFL = ItF; k2=2; k4=4; j1=1; j4=4; m2=2; m3=3;
case 2
 pos_ = and(nflp,ItA); ItFL = ItF; k2=4; k4=2; j1=4; j4=1; m2=3; m3=2;
case 3
 pos_ = and(flp,ItB); ItFL = ItL; k2=2; k4=4; j1=1; j4=4; m2=2; m3=3;
case 4
 pos_ = and(nflp,ItB); ItFL = ItL; k2=4; k4=2; j1=4; j4=1; m2=3; m3=2;
end;
if not(any(pos_))
continue;
end;
 pos = repmat(pos_,2,1); 
 ItF_ = and(pos,ItFL); ItL_ = and(pos,not(ItFL));
 newtri3_ = [[fac(ItL_,2) fac(ItF_,[j4 j1 3])]; ...
             [fac(ItL_,2) fac(ItF_,[j1 j4 2])]; ...
             [fac(ItF_,2) fac(ItL_,[m3 m2 1])]; ...
             [fac(ItF_,2) fac(ItL_,[m2 m3 4])]; ...
             [fac(ItL_,[2 4]) fac(ItF_,[k2 k4])]];
 %[tmp1,vol1] =elem_inv(newtri3_,xy); vol1=reshape(vol1,sum(ItF_),5); tmp1=reshape(tmp1,sum(ItF_),5);
 %[tmp,vol] =elem_inv(tri(I,:),xy); vol=vol(nEA); %plot(tmp1)
 %[sum(vol(pos(ItF))) sum(vol1(:))]
 %Call_ = Call(pos(ItF)); unique(Call_(not(any(tmp1,2))))
 newtri3 = [newtri3; newtri3_];
 oldtri = [oldtri; repmat(fI3n(pos_),5,1)];
 q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
end; %fori
end; %ifany
end; %ifany nEA
%ItFL = [thmsk(Call);not(thmsk(Call))];
%for i=1:2
%if i == 1
 %pos = It1; m1=1; m4=4; j2=2; j3=3; k2=2; k4=4; 
%else
 %pos = It2; m1=2; m4=3; j2=4; j3=1; k4=4; k2=2;
%end;
%if not(any(pos))
%continue;
%end;
 %pos = repmat(pos,2,1);
 %ItF_ = and(pos,ItFL); ItL_ = and(pos,not(ItFL));
 %TriT1 = [[fac(ItF_,2) fac(ItL_,[m4 m1 j3])]; ...
         %[fac(ItL_,2) fac(ItF_,[m4 m1 j3])]; ...
         %[fac(ItF_,2) fac(ItL_,[m1 m4 j2])]; ...
         %[fac(ItL_,2) fac(ItF_,[m1 m4 j2])]; ...
         %[fac(ItL_,[2 4]) fac(ItF_,[k2 k4])]];
 %if i==1
 %N3_ = sum(ItF_);
 %TriT2 = [[fac(ItF_,2) fac(ItL_,[m1 3 m4])]; ...
         %[fac(ItL_,2) fac(ItF_,[m1 3 m4])]; ...
         %[fac(ItF_,4) fac(ItL_,[m4 2 m1])]; ...
         %[fac(ItL_,4) fac(ItF_,[m4 2 m1])]];
 %qual1 = elem_qual(TriT1,xy,Nmetric,options);
 %qual2 = elem_qual(TriT2,xy,Nmetric,options);
 %qual1m = min(reshape(qual1,N3_,5),[],2);
 %qual2m = min(reshape(qual2,N3_,4),[],2);
 %option1 = qual1m > qual2m; option2 = not(option1);
 %option1 = repmat(option1,5,1);
 %option2 = repmat(option2,4,1);
 %newtri3_ = [TriT1(option1,:); TriT2(option2,:)];
 %q3 = [q3; qual1(option1,:); qual2(option2,:)];
 %else
 	%newtri3_ = TriT1;
 	%q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
 %end;
 %%[tmp1,vol1] =elem_inv(newtri3_,xy); vol1=reshape(vol1,sum(ItF_),4); tmp1=reshape(tmp1,sum(ItF_),4);
 %%[tmp,vol] =elem_inv(tri(I,:),xy); vol=vol(nEA); %plot(tmp1)
 %%[sum(vol(pos(ItF))) sum(vol1(:))]
  %%[vol(pos(ItF)) sum(vol1,2)]
 %% [Call(pos(ItF)) any(tmp1,2)]
 %%Call_ = Call(pos(ItF)); unique(Call_(not(any(tmp1,2))))
 %newtri3 = [newtri3; newtri3_];
%end; %fori
%
%%no symmetry
%ItA = and(qsm1(ItF),nqsm1(ItL));
%ItB = and(nqsm1(ItF),qsm1(ItL));
%if any([ItA(:);ItB(:)])
%for i=1:2
%if i == 1
 %pos = ItA; k2=2; k4=4; j1=1; j4=4; m2=2; m3=3;
%else
 %pos = ItB; ItFL = not(ItFL); k2=2; k4=4; j1=1; j4=4; m2=2; m3=3;
%end;
%if not(any(pos))
%continue;
%end;
 %pos = repmat(pos,2,1); 
 %ItF_ = and(pos,ItFL); ItL_ = and(pos,not(ItFL));
 %newtri3_ = [[fac(ItL_,2) fac(ItF_,[j4 j1 3])]; ...
             %[fac(ItL_,2) fac(ItF_,[j1 j4 2])]; ...
             %[fac(ItF_,2) fac(ItL_,[m3 m2 1])]; ...
             %[fac(ItF_,2) fac(ItL_,[m2 m3 4])]; ...
             %[fac(ItL_,[2 4]) fac(ItF_,[k2 k4])]];
 %%[tmp1,vol1] =elem_inv(newtri3_,xy); vol1=reshape(vol1,sum(ItF_),5); tmp1=reshape(tmp1,sum(ItF_),5);
 %%[tmp,vol] =elem_inv(tri(I,:),xy); vol=vol(nEA); %plot(tmp1)
 %%[sum(vol(pos(ItF))) sum(vol1(:))]
 %%Call_ = Call(pos(ItF)); unique(Call_(not(any(tmp1,2))))
 %newtri3 = [newtri3; newtri3_];
 %q3 = [q3; elem_qual(newtri3_,xy,Nmetric,options)];
%end; %fori
%end; %ifany
%end; %ifany nEA
%[tmp,vol] =elem_inv(tri(I,:),xy);
%[tmp1,vol1] =elem_inv(newtri3,xy);
%[sum(vol) sum(vol1)]


function [newtri4,q4,xy,Nmetric,oldtri] = mesh4long(deltri4_,xy,edga,tri2edga4,fac2,fac2q,delfacc2,tri2fac4,fac2tri2,fI4,Nmetric,options)
%deltri4_ = deltri_(I,:);
%tri2edga4 = tri2edga(I,:);
tzr = max(edga(:)); 
Ca1  = [1 3 5 6 9 11]'; %find([true false true false true true false false true false true false]) (gen_books3D)
EA = any([and(deltri4_(:,1)==tzr,deltri4_(:,2)==tzr) ...
          and(deltri4_(:,3)==tzr,deltri4_(:,4)==tzr) ...
          and(deltri4_(:,5)==tzr,deltri4_(:,6)==tzr)],2); nEA = not(EA);
          
%CASE 1, one isolated node
if any(nEA)
deltri4n_ = deltri4_(nEA,:);
tri2edga4n = tri2edga4(nEA,:); 
fI4n = fI4(nEA); N4 = size(deltri4n_,1);

nR = (1:N4)';
CE = 2*all(deltri4n_(:,[1 2]) ~= tzr,2) ...
   + 4*all(deltri4n_(:,[3 4]) ~= tzr,2) ...
   + 6*all(deltri4n_(:,[5 6]) ~= tzr,2); %split edge pair
edgaC = [CE CE-1];
crnrs = [edga(tri2edga4n(nR+(Ca1(edgaC(:,1))-1)*N4),:) ...
         edga(tri2edga4n(nR+(Ca1(edgaC(:,2))-1)*N4),:)];
Cs = repmat((1:6)',1,N4); Cs = reshape(Cs(deltri4n_'==tzr),2,N4)';
thmsk = [0 0 1 4 1 3; ...
	0 0 3 2 4 2; ...
	0 0 0 0 2 4; ...
	0 0 0 0 3 1];
Cs_ = thmsk(Cs(:,1)+(Cs(:,2)-1)*4);
%ndlock = crnrs(nR+(Cs_-1)*N4);
thmsk = [2 1 4 3]'; 
loneCR = crnrs(nR+(thmsk(Cs_)-1)*N4);
thmsk = [2 2 1 1]'; edgaC2 = edgaC(nR+(thmsk(Cs_)-1)*N4);
O5n  = deltri4n_(nR+(edgaC2-1)*N4);
ItF = [true(N4,1); false(N4,1)]; ItL = not(ItF);

tri2fac4t = tri2fac4(nEA,:)'; fact = delfacc2(tri2fac4t);
fact_ = tri2fac4t; fact_ = reshape(tri2fac4t(fact~=0),2,N4)';
trinEA_ = reshape(fac2tri2(fact_),N4,2); %reshape  relevant for N4==1
fact = reshape(fact(fact~=0),2,N4)'; qSm1 = fac2q(fact(:)); facO = fac2(fact(:),1:4); ndlock = fac2(fact(:,1),5);
flp_ = trinEA_ == repmat(fI4n,1,2);
facO(flp_(:),:) = facO(flp_(:),[2 1 4 3]); qSm1(flp_(:)) = not(qSm1(flp_(:)));

%enforce facO(ItF,1) == fac(ItL,1)
flp = facO(ItL,1) ~= loneCR; ItF_ = [flp; false(N4,1)]; ItL_ = [false(N4,1); flp];
tmp = facO(ItF_,:); facO(ItF_,:) = facO(ItL_,:); facO(ItL_,:) = tmp;
tmp = qSm1(ItF_); qSm1(ItF_) = qSm1(ItL_); qSm1(ItL_) = tmp;
facO(ItF,:) = facO(ItF,[2 1 4 3]);
qSm1(ItF) = not(qSm1(ItF));
nqSm1 = not(qSm1); facO = facO(:,[1 3 4 2]); nqSm1 = qSm1; qSm1 = not(nqSm1);%could be avoided

%%%brute find
%%CE1 = CE-2; CE1(CE1==0) = 6;
%%CE2 = CE+2; CE2(CE2==8) = 2;
%%nedga1 = [deltri4n_(nR+(CE1-2)*N4) deltri4n_(nR+(CE1-1)*N4)]; nedga1 = nedga1(:,1)~=tzr; 
%%nedga2 = [deltri4n_(nR+(CE2-2)*N4) deltri4n_(nR+(CE2-1)*N4)]; nedga2 = nedga2(:,1)~=tzr;
%%CE1(nedga1) = CE1(nedga1)-1; nedga1 = tri2edga4n(nR+(Ca1(CE1)-1)*N4);
%%CE2(nedga2) = CE2(nedga2)-1; nedga2 = tri2edga4n(nR+(Ca1(CE2)-1)*N4);
%%ndlock_ = sort([edga(nedga1,:) edga(nedga2,:)],2)';
%%ndlock_ = ndlock_([false(1,size(nedga1,1)); ndlock_(1:3,:)==lonen_(2:4,:)]);
%%ndlock==ndlock_
%
%thmsk = [1 1 2 2]'; edgaC1 = edgaC(nR+(thmsk(Cs_)-1)*N4);
%ndlone = deltri4n_(nR+(edgaC1-1)*N4);
%Cs2 = deltri4n_~=tzr; Cs2([nR; nR]+(edgaC(:)-1)*N4) = 0;
%thmsk = [[3 4 1 2]; ...
	%[4 3 2 1]]'; Cs2_ = [thmsk(Cs_,1); thmsk(Cs_,2)];
%deltri4n_t = deltri4n_'; sides = [reshape(deltri4n_t(Cs2'),2,N4)' reshape(crnrs([nR; nR]+(Cs2_-1)*N4),N4,2)];
%flpI = and(or(Cs(:,1)==1,Cs(:,1)==2),or(Cs(:,2)==5,Cs(:,2)==6));
%sides(flpI,:) = sides(flpI,[2 1 3 4]);
%flpI = or(Cs(:,1)==3,Cs(:,1)==4);
%sides(flpI,:) = sides(flpI,[1 2 4 3]);
%
%%[tmp1,vol1] = elem_inv(newtri4,xy); %[tmp1 Cs]
    %%[tmp1,vol] = elem_inv(tri(I,:),xy);
    %%[sum(vol(nEA)) sum(vol1)*8]
%%[tmp1,area1] = elem_inv(sides(:,[1 2 3]),xy);
%%[tmp2,area2] = elem_inv(sides(:,[2 3 4]),xy); area12 = area1+area2;
%%[tmp,area] = elem_inv([nlone sides(:,[3 4])],xy);
%%area./area12*3/4
%%ijk=2; l1 = [ndlock(ijk) sides(ijk,1) O5n(ijk) sides(ijk,2) ndlock(ijk)]; 
%%l2a = [ndlone(ijk) sides(ijk,1) ndlock(ijk) loneCR(ijk) sides(ijk,4)];
%%l2b = [ndlone(ijk) sides(ijk,2) ndlock(ijk) loneCR(ijk) sides(ijk,3)]; figure(); plot3(xy(l1,1),xy(l1,2),xy(l1,3),'-b', xy(l2a,1),xy(l2a,2),xy(l2a,3),'-m', xy(l2b,1),xy(l2b,2),xy(l2b,3),'-r');
%
%facO = [[loneCR;loneCR] [ndlone;ndlone] [sides(:,[1 4]); sides(:,[2 3])]];
newtri4 = [facO(ItF,[3 2]) facO(ItL,3) ndlock];
q4 = elem_qual(newtri4,xy,Nmetric,options);
oldtri = fI4n;
%
%%qualS1 = elem_qual([facO(:,[2 1 3]); ...
		  %%facO(:,[4 1 3])],xy,Nmetric,options);
%%qualS2 = elem_qual([facO(:,[1 2 4]); ...
		  %%facO(:,[3 2 4])],xy,Nmetric,options);
%%qualS1m = min(reshape(qualS1,N4*2,2),[],2);
%%qualS2m = min(reshape(qualS2,N4*2,2),[],2);
%%qSm1 = qualS1m > qualS2m; nqSm1 = not(qSm1);
%[qSm1,nqSm1] = surfqual(facO,[[2 1 3];[4 1 3]],[[1 2 4]; [3 2 4]],xy,Nmetric,options); 

%all to facO(:,1) (25% of cases)
It = all(reshape(qSm1,N4,2),2);
if any(It)
 ItF_ = [It; false(N4,1)]; ItL_ = [false(N4,1); It];
 O5n_ = O5n(It);
 newtri4_ = [[facO(ItF_,[1 3]) facO(ItL_,[2 3])]; ...
 	    [facO(ItF_,[1 3]) facO(ItL_,3) O5n_]; ...
 	    [facO(ItF_,[1 4 3]) O5n_]; ...
 	    [facO(ItL_,[1 3 4]) O5n_]];
 %[tmp1,vol1] = elem_inv(newtri4_,xy); [tmp1 repmat(Cs(It,1),4,1)]
 %[tmp,vol] = elem_inv(tri(I,:),xy); vol=vol(nEA);
 %[sum(vol(It)) sum(vol1)*8/7]
 newtri4 = [newtri4; newtri4_];
 q4 = [q4; elem_qual(newtri4_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(fI4n(It),4,1)];
end; %ifany

%all to facO(:,2) (25% of cases)
It = all(reshape(nqSm1,N4,2),2);
if any(It)
 ItF_ = [It; false(N4,1)]; ItL_ = [false(N4,1); It];
 O5n_ = O5n(It);
 newtri4_ = [[facO(ItF_,[2 4 3]) O5n_]; ...
 	    [facO(ItL_,[2 3 4]) O5n_]; ...
 	    [facO(ItF_,[1 4 2]) O5n_]; ...
 	    [facO(ItL_,[1 2 4]) O5n_]; ...
 	    [facO(ItF_,[2 3]) facO(ItL_,3) O5n_]];
 %[tmp1,vol1] = elem_inv(newtri4_,xy); [tmp1 repmat(Cs(It,1),5,1)]
 %[tmp1,vol] = elem_inv(tri(I,:),xy); vol=vol(nEA);
 %[sum(vol(It)) sum(vol1)*8/7]
 newtri4 = [newtri4; newtri4_];
 q4 = [q4; elem_qual(newtri4_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(fI4n(It),5,1)];
end; %ifany

%mixed cases (last 50%)
It=or(and(qSm1(ItF),nqSm1(ItL)),and(nqSm1(ItF),qSm1(ItL)));
if any(It)
It1=and(qSm1(ItF),nqSm1(ItL)); %first to facO(:,1)
It2=and(nqSm1(ItF),qSm1(ItL)); %first to facO(:,2)
for i=1:2
 if i==1
  pos_ = It1; ItFL = ItF; j1=1; j2=2; j3=3; j4=4;
 else
  pos_ = It2; ItFL = ItL; j1=4; j2=3; j3=2; j4=1;
 end;
 if not(any(pos_))
 	continue;
 end;
 O5n_ = O5n(pos_);
 pos = repmat(pos_,2,1);
 ItF_ = and(pos,ItFL); ItL_ = and(pos,not(ItFL));
 newtri4_ = [[facO(ItF_,[1 j3 j2]) O5n_]; ...
 	    [facO(ItF_,[j4 3 j1]) O5n_]; ...
 	    [facO(ItL_,[j3 4 j2]) O5n_]; ...
 	    [facO(ItL_,[j4 j1 2]) O5n_]; ...
 	    [facO(ItF_,[j2 j3]) facO(ItL_,3) O5n_]];
 %[tmp1,vol1] = elem_inv(newtri4_,xy); [tmp1 repmat(Cs(pos(ItF),1),5,1)]
 %[tmp1,vol] = elem_inv(tri(I,:),xy); vol=vol(nEA);
 %[sum(vol(pos(ItF))) sum(vol1)*8/7]
 newtri4 = [newtri4; newtri4_];
 q4 = [q4; elem_qual(newtri4_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(fI4n(pos_),5,1)];
end; %fori
end; %ifany
%[tmp1,vol1] = elem_inv(newtri4,xy);
%[tmp1,vol] = elem_inv(tri(I,:),xy); vol=vol(nEA);
%[sum(vol) sum(vol1)]
else %ifany EA
 newtri4 = [];  q4 = []; oldtri = [];
end; %ifany EA

%CASE 2, no isolated nodes
if any(EA)
 deltri4n_ = deltri4_(EA,:);
 tri2edga4n = tri2edga4(EA,:); 
 fI4n = fI4(EA); N4 = size(deltri4n_,1);
 ItF = [true(N4,1); false(N4,1)]; ItL = not(ItF);
 %triE = tri(I,:);
 %crnrs = triE(EA,:);
 [CE,Rn] = find(deltri4n_(:,[1 3 5])'==tzr);
 
 %edgaC = [2*CE-1 2*CE];
 %crnrs2 = [edga(tri2edga4n(Rn+(Ca1(edgaC(:,1))-1)*N4),:) ...
 %          edga(tri2edga4n(Rn+(Ca1(edgaC(:,2))-1)*N4),:)];
 deltri4n_t = deltri4n_'; 
 CEp = CE+1; CEp(CEp==4) = 1;
 %CEm = CE-1; CEm(CEm==0) = 3;
 edgaC = [2*CEp-1 2*CEp];
 crnrs3 = [edga(tri2edga4n(Rn+(Ca1(edgaC(:,1))-1)*N4),:) ...
           edga(tri2edga4n(Rn+(Ca1(edgaC(:,2))-1)*N4),:)];
 crnrs3(CE~=3,:) = crnrs3(CE~=3,[2 1 3 4]);

 tri2fac4t = tri2fac4(EA,:)'; fact = delfacc2(tri2fac4t);
fact_ = tri2fac4t; fact_ = reshape(tri2fac4t(fact~=0),4,N4)';
%trinEA_ = fac2tri2(fact_);
factt = reshape(fact(fact~=0),4,N4); fact = factt'; ndcr = reshape(fac2(fact(:),5),N4,4)';
tri2fac4t = repmat(tri2fac4t,1,2); factt = repmat(factt,1,2); fI4_ = repmat(fI4n,2,1);
iA = [repmat(crnrs3(:,1)',4,1) repmat(crnrs3(:,4)',4,1)] == repmat(ndcr,1,2); triA = fac2tri2(tri2fac4t(iA)); factA = factt(iA); facA = fac2(factA,1:4); qAm1 = fac2q(factA);
iB = [repmat(crnrs3(:,2)',4,1) repmat(crnrs3(:,3)',4,1)] == repmat(ndcr,1,2); triB = fac2tri2(tri2fac4t(iB)); factB = factt(iB); facB = fac2(factB,1:4); qBm1 = fac2q(factB);
flpA = triA == fI4_; facA(flpA,:) = facA(flpA,[2 1 4 3]); qAm1(flpA) = not(qAm1(flpA));
flpB = triB == fI4_; facB(flpB,:) = facB(flpB,[2 1 4 3]); qBm1(flpB) = not(qBm1(flpB));

%enforce facA(ItF,[1 2]) == facA(ItL,[1 2])
facA(ItF,:) = facA(ItF,[2 1 4 3]); qAm1(ItF) = not(qAm1(ItF));
facB(ItF,:) = facB(ItF,[2 1 4 3]); qBm1(ItF) = not(qBm1(ItF));
nqAm1 = qAm1; facA = facA(:,[1 2 4 3]); qAm1 = not(nqAm1);
nqBm1 = qBm1; facB = facB(:,[1 2 4 3]); qBm1 = not(nqBm1);%could be avoided
 
  
 %sides = reshape(deltri4n_t(deltri4n_t~=tzr),4,N4)';
 %sides(CE==1,:) = sides(CE==1,[1 2 4 3]);
 %sides(CE==2,:) = sides(CE==2,[2 1 3 4]);
 %sides(CE~=2,:) = sides(CE~=2,[3 4 1 2]);
 %%[tmp,vol] = elem_inv(crnrs3,xy); [tmp CE]
 %%[tmp,area] = elem_inv(crnrs3(:,[2 3 4]),xy);
 %%ijk = 4;
 %%[tmp,area1] = elem_inv([crnrs3(:,ijk) sides(:,[1 4])],xy); sum(abs(area./area1/4-1)<1e-6) %2
 %%[tmp,area1] = elem_inv([crnrs3(:,ijk) sides(:,[2 3])],xy); sum(abs(area./area1/4-1)<1e-6)
 %%[tmp,area1] = elem_inv([crnrs3(:,ijk) sides(:,[1 3])],xy); sum(abs(area./area1/4-1)<1e-6)
 %%[tmp,area1] = elem_inv([crnrs3(:,ijk) sides(:,[2 4])],xy); sum(abs(area./area1/4-1)<1e-6)
 %facA = [[crnrs3(:,[2 3]); crnrs3(:,[2 3])] [sides(:,[1 3]); sides(:,[4 2])]];
 %facB = [[crnrs3(:,[1 4]); crnrs3(:,[1 4])] [sides(:,[2 3]); sides(:,[4 1])]];
%% figure(); ijk=1; ijk2 = ijk+N4; plot3(xy(facA(ijk,:),1),xy(facA(ijk,:),2),xy(facA(ijk,:),3),'-k',xy(facB(ijk,:),1),xy(facB(ijk,:),2),xy(facB(ijk,:),3),'-r', xy(facA(ijk2,:),1),xy(facA(ijk2,:),2),xy(facA(ijk2,:),3),'-b',xy(facB(ijk2,:),1),xy(facB(ijk2,:),2),xy(facB(ijk2,:),3),'-m')
 %%qualA1 = elem_qual([facA(:,[2 1 3]); ...
 		   %%facA(:,[4 1 3])],xy,Nmetric,options);
 %%qualA2 = elem_qual([facA(:,[1 2 4]); ...
 		   %%facA(:,[3 2 4])],xy,Nmetric,options);
 %%qualB1 = elem_qual([facB(:,[2 1 3]); ...
 		   %%facB(:,[4 1 3])],xy,Nmetric,options);
 %%qualB2 = elem_qual([facB(:,[1 2 4]); ...
 		   %%facB(:,[3 2 4])],xy,Nmetric,options);
%%qualA1m = min(reshape(qualA1,N4*2,2),[],2);
%%qualA2m = min(reshape(qualA2,N4*2,2),[],2);
%%qualB1m = min(reshape(qualB1,N4*2,2),[],2);
%%qualB2m = min(reshape(qualB2,N4*2,2),[],2);
%%qAm1 = qualA1m > qualA2m; nqAm1 = not(qAm1);
%%qBm1 = qualB1m > qualB2m; nqBm1 = not(qBm1);
%[qAm1,nqAm1] = surfqual(facA,[[2 1 3];[4 1 3]],[[1 2 4]; [3 2 4]],xy,Nmetric,options); 
%[qBm1,nqBm1] = surfqual(facB,[[2 1 3];[4 1 3]],[[1 2 4]; [3 2 4]],xy,Nmetric,options); 

%all to something (25 % of cases)
It1 = all(reshape(nqAm1,N4,2),2);
It2 = all(reshape(qAm1,N4,2),2);
It3 = all(reshape(nqBm1,N4,2),2);
It4 = all(reshape(qBm1,N4,2),2);
for i=1:4
 switch i
 case 1 %THIS ALWAYS TRIGGERS CONTINUE!?
  pos_=and(It1,It3); k1=2; k2=1; j4=4; %to facA(:,2) and facB(:,2)
  ItFL = ItF;
 case 2
  pos_=and(It1,It4); k1=1; k2=2; j4=4; %to facA(:,2) and facB(:,1)
  ItFL = ItL;
 case 3
  pos_=and(It2,It3); k1=2; k2=1; j4=3; %to facA(:,1) and facB(:,2)
  ItFL = ItF;
 otherwise
  pos_=and(It2,It4); k1=1; k2=2; j4=3; %to facA(:,1) and facB(:,1)
  ItFL = ItL;
 end;
 if not(any(pos_))
  continue;
 end;
 N4_ = sum(pos_); j2=j4-2; k4=k1+2;
 ItF_ = [pos_; false(N4,1)]; ItL_ = [false(N4,1); pos_];
 pos = repmat(pos_,2,1);
 ItF2_ = and(pos,ItFL); ItL2_ = and(pos,not(ItFL));
 TriT1 = [[facA(ItF_,[j2 3 4]) facA(ItL_,3)]; ...
	 [facA(ItL_,[j2 4 3]) facA(ItF_,4)]; ...
	 [facB(ItF_,[k1 3 4]) facB(ItL_,3)]; ...
	 [facB(ItL_,[k1 4 3]) facB(ItF_,4)]]; %4F to 3L
 TriT2 = [[facA(ItF_,[j2 3 4]) facA(ItL_,4)]; ...
	 [facA(ItL_,[j2 4 3]) facA(ItF_,3)]; ...
	 [facB(ItF_,[k1 3 4]) facB(ItL_,4)]; ...
	 [facB(ItL_,[k1 4 3]) facB(ItF_,3)]]; %A3F to A4L or (B4L to B3F)
 %[tmp1,vol] = elem_inv(TriT1,xy);
 %[tmp2,vol] = elem_inv(TriT2,xy); 
 qual1 = elem_qual(TriT1,xy,Nmetric,options);
 qual2 = elem_qual(TriT2,xy,Nmetric,options);
 qual1m = min(reshape(qual1,N4_,4),[],2);
 qual2m = min(reshape(qual2,N4_,4),[],2);
 option1 = qual1m > qual2m; option2 = not(option1); 

 option1 = repmat(option1,4,1); 
 option2 = repmat(option2,4,1); 
 %newtri4_ = zeros(N4_,4); newtri4_(option1,:) = TriT1(option1,:); newtri4_(option2,:) = TriT2(option2,:);
 newtri4_ = [TriT1(option1,:); TriT2(option2,:)];
 q4_ = [qual1(option1); qual2(option2)];
 oldtri_ = repmat(fI4n(pos_),4,1); oldtri_ = [oldtri_(option1); oldtri_(option2)];
 newtri4__ = [[facA(ItF_,[1 2 j4]) facA(ItL_,j4)]; ...
	     [facB(ItF2_,[k2 k1 k4]) facB(ItL2_,k4)]];
 %[tmp1,vol1] = elem_inv(newtri4__,xy); [tmp1 repmat(CE(ItF_),2,1)]
 newtri4 = [newtri4; newtri4_; newtri4__];
 q4 = [q4; q4_; elem_qual(newtri4__,xy,Nmetric,options)];
 oldtri = [oldtri; oldtri_; repmat(fI4n(pos_),2,1)];

 %[tmp1,vol1] = elem_inv(newtri4_,xy); 
 %[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA);
 %[i sum(vol1) sum(vol(pos(ItF)))]
end; %fori

%mixed cases (next 62.5 % of cases)
ItA = and(nqAm1(ItF),qAm1(ItL)); %A4F to A2, A3L to A1
ItB = and(qAm1(ItF),nqAm1(ItL)); %A3F to A1, A4L to A2
ItC = and(nqBm1(ItF),qBm1(ItL)); %B4F to B2, B3L to B1
ItD = and(qBm1(ItF),nqBm1(ItL)); %B3F to B1, B4L to B2
for i=1:10
 switch i
 case 1
 pos_=and(It1,ItC); j4=4; ItFL = ItF; m3=3; m4=4; %4F to 3L
 case 2
 pos_=and(It1,ItD); j4=4; ItFL = ItL; m3=4; m4=3; %A3F=B4L to A4L=B3F
 case 3
 pos_=and(It2,ItC); j4=3; ItFL = ItF; m3=3; m4=4; %4F to 3L
 case 4
 pos_=and(It2,ItD); j4=3; ItFL = ItL; m3=4; m4=3; %A3F=B4L to A4L=B3F
 case 5
 pos_=and(It3,ItA); k1=2; ItFL = ItF; m3=3; m4=4; %4F to 3L
 case 6
 pos_=and(It3,ItB); k1=2; ItFL = ItL; m3=4; m4=3;
 case 7
 pos_=and(It4,ItA); k1=1; ItFL = ItF; m3=3; m4=4; %4F to 3L
 case 8
 pos_=and(It4,ItB); k1=1; ItFL = ItL; m3=4; m4=3;
 case 9 %supermix cases
 pos_=and(ItA,ItC); ItFL = ItF; m3=3; m4=4; %4F to 3L
 otherwise
 pos_=and(ItB,ItD); ItFL = ItL; m3=4; m4=3; %A3F=B4L to A4L=B3F
 end; %[i sum(pos)]
 if not(any(pos_))
 	continue;
 end;
 ItF_ = [pos_; false(N4,1)]; ItL_ = [false(N4,1); pos_]; 
 pos = repmat(pos_,2,1); m1=m3-2; m2=m4-2;
 ItF2_ = and(pos,ItFL); ItL2_ = and(pos,not(ItFL));
 if i==1 || i==2 || i==3 || i==4
  j2=j4-2; 
  newtri4_ = [[facA(ItF_,[1 2 j4]) facA(ItL_,j4)]; ...
              [facA(ItF_,[j2 3 4]) facA(ItL_,m3)]; ...
	     [facA(ItL_,[j2 4 3]) facA(ItF_,m4)]; ...
	     [facB(ItF2_,[2 m3 m4]) facB(ItL2_,3)]; ...
	     [facB(ItF2_,[m2 4 m1]) facB(ItL2_,3)]; ...
	     [facB(ItL2_,[1 m4 m3]) facB(ItF2_,4)]];
 elseif i==5 || i==6 || i==7 || i==8
  k4=k1+2; 
  newtri4_ = [[facB(ItF_,[1 2 k4]) facB(ItL_,k4)]; ...
	     [facB(ItF_,[k1 3 4]) facB(ItL_,m3)]; ...
	     [facB(ItL_,[k1 4 3]) facB(ItF_,m4)]; ...
	     [facA(ItF2_,[2 m3 m4]) facA(ItL2_,3)]; ...
	     [facA(ItF2_,[m2 4 m1]) facA(ItL2_,3)]; ...
	     [facA(ItL2_,[1 m4 m3]) facA(ItF2_,4)]];
 else
  newtri4_ = [[facA(ItF2_,[2 m3 m4]) facA(ItL2_,3)]; ...
	     [facA(ItF2_,[m2 4 m1]) facA(ItL2_,3)]; ...
	     [facA(ItL2_,[1 m4 m3]) facA(ItF2_,4)]; ...
	     [facB(ItF2_,[2 m3 m4]) facB(ItL2_,3)]; ...
	     [facB(ItF2_,[m2 4 m1]) facB(ItL2_,3)]; ...
	     [facB(ItL2_,[1 m4 m3]) facB(ItF2_,4)]];
 end;
 newtri4 = [newtri4; newtri4_];
 q4 = [q4; elem_qual(newtri4_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(fI4n(pos_),6,1)];

 %[tmp1,vol1] = elem_inv(newtri4_,xy);
 %[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA);
 %[i sum(vol1) sum(vol(pos(ItF)))]
end; %fori

%impossible mixed cases (12.5 %)
It1=and(ItA,ItD);
It2=and(ItB,ItC); 
%if any([It1;It2])
	%save for_debug3D.mat; error('YE, found input to test untested code'); %error('just stop 3D');
%end;
if any([It1;It2])
for i=1:2
 if i==1
 pos_ = It1; ItFL = ItF; j3=3; j4=4; k1=1; k4=4; k2=2; k3=3;
 else
 pos_ = It2; ItFL = ItL; j3=4; j4=3; k1=4; k4=1; k2=3; k3=2;
 end;
 if not(any(pos_))
 	continue;
 end;
 In = [facA(pos_,1) facA(pos_,2) facB(pos_,1) facB(pos_,2)];
 nxy = 0.25*(xy(In(:,1),:)+xy(In(:,2),:)+xy(In(:,3),:)+xy(In(:,4),:));
 nNmetric =  metric_avg(In,Nmetric,options);
 In = (1:size(nxy,1))'+size(xy,1);
 pos = repmat(pos_,2,1); j2=j4-2; j1=j3-2; 
 ItF2_ = and(pos,ItFL); ItL2_ = and(pos,not(ItFL));
 newtri4_ = [[In facA(ItF2_,[2 j4 j3])]; ...
	    [In facA(ItF2_,[j2 j1 4])]; ...
 	    [In facA(ItL2_,[1 j3 j4])]; ...
	    [In facA(ItL2_,[j1 j2 3])]; ...
	    [In facA(ItF2_,[k4 k1]) facA(ItL2_,4)]; ...
	    [In facA(ItF2_,[k2 k3]) facA(ItL2_,3)]; ...
	    [In facB(ItL2_,[2 j3 j4])]; ...
	    [In facB(ItL2_,[j1 j2 4])]; ...
 	    [In facB(ItF2_,[1 j4 j3])]; ...
	    [In facB(ItF2_,[j2 j1 3])]; ...
	    [In facB(ItL2_,[k1 k4]) facB(ItF2_,4)]; ...
	    [In facB(ItL2_,[k3 k2]) facB(ItF2_,3)]];
 xy_ = [xy; nxy];
 %[tmp1,vol1] = elem_inv(newtri4_,xy_);  %any(tmp1)
 %[tmp,vol] = elem_inv(tri(I,:),xy); vol = vol(EA);
 %[i sum(vol1) sum(vol(pos(ItF)))]
 xy = [xy; nxy];
 Nmetric = [Nmetric; nNmetric];
 newtri4 = [newtri4; newtri4_];
 q4 = [q4; elem_qual(newtri4_,xy,Nmetric,options)];
 oldtri = [oldtri; repmat(fI4n(pos_),12,1)];

 %figure(); N4_ = sum(ItF2_); ijk=1; ptri=newtri4_(ijk:N4_:end,:); 
 %plot3(reshape(xy(ptri(:,[1 2]),1),size(ptri,1),2)', reshape(xy(ptri(:,[1 2]),2),size(ptri,1),2)', reshape(xy(ptri(:,[1 2]),3),size(ptri,1),2)','-b', reshape(xy(ptri(:,[1 3]),1),size(ptri,1),2)', reshape(xy(ptri(:,[1 3]),2),size(ptri,1),2)', reshape(xy(ptri(:,[1 3]),3),size(ptri,1),2)','-b', reshape(xy(ptri(:,[1 4]),1),size(ptri,1),2)', reshape(xy(ptri(:,[1 4]),2),size(ptri,1),2)', reshape(xy(ptri(:,[1 4]),3),size(ptri,1),2)','-b', reshape(xy(ptri(:,[2 3]),1),size(ptri,1),2)', reshape(xy(ptri(:,[2 3]),2),size(ptri,1),2)', reshape(xy(ptri(:,[2 3]),3),size(ptri,1),2)','-b' ,reshape(xy(ptri(:,[2 4]),1),size(ptri,1),2)', reshape(xy(ptri(:,[2 4]),2),size(ptri,1),2)', reshape(xy(ptri(:,[2 4]),3),size(ptri,1),2)','-b', reshape(xy(ptri(:,[3 4]),1),size(ptri,1),2)', reshape(xy(ptri(:,[3 4]),2),size(ptri,1),2)', reshape(xy(ptri(:,[3 4]),3),size(ptri,1),2)','-b');
end; %fori
end; %anyif
end; %ifany EA
%[tmp4,vol4] = elem_inv(newtri4,xy); 
%[tmp,vol] = elem_inv(tri(I,:),xy);
%[sum(vol) sum(vol4)]



function [newtri5,q5,oldtri] = mesh5long(deltri5_,xy,edga,tri2edga5,fac2,fac2q,delfacc2,tri2fac5,fac2tri2,fI5,Nmetric,options)
    %tri2edga5 = tri2edga(I,:);
    %deltri5_ = deltri_(I,:);
    N5 = size(deltri5_,1);
    tzr = max(edga(:));
    [C,R] = find(deltri5_'==tzr);
    Copp = [2 1 4 3 6 5]';
    ndshr = deltri5_(R+(Copp(C)-1)*N5);
    Ca1  = [1 3 5 6 9 11]'; %find([true false true false true true false false true false true false]) (gen_books3D)
    ndfreeA = edga(tri2edga5(R+(Ca1(C)-1)*N5),:);
    ndlockA = edga(tri2edga5(R+(Ca1(Copp(C))-1)*N5),:);
    Cp2 = C  +2; Cp2(Cp2==7) = 1; Cp2(Cp2==8) = 2;
    Cp4 = Cp2+2; Cp4(Cp4==7) = 1; Cp4(Cp4==8) = 2;
    others = [deltri5_(R+(Cp2-1)*N5) deltri5_(R+(Copp(Cp2)-1)*N5) ... 
    deltri5_(R+(Cp4-1)*N5) deltri5_(R+(Copp(Cp4)-1)*N5)];
    It2_ = any([C==2 C==4 C==6],2); 
    others(It2_,[3 4]) = others(It2_,[4 3]);
    It3_ = C==1;
    others(It3_,:) = others(It3_,[2 1 4 3]);
    %It4_ = C==5; 
    %others(It4_,:) = others(It4_,[2 1 4 3]);
    newtri5 = zeros(N5*2,4);
    newtri5(:,1) = [ndshr;ndshr];
    newtri5(:,2) = ndlockA(:);
    newtri5(:,[3 4]) = [others(:,[4 2]);others(:,[3 1])];
    q5 = elem_qual(newtri5,xy,Nmetric,options);
    oldtri = repmat(fI5,2,1);
    %[tmp5,vol5] = elem_inv(newtri5,xy);  [tmp,vol2] = elem_inv(tri(I,:),xy); [sum(vol2)-sum(vol5)*4 all(vol5>0)]
    %unique(C(abs(abs(vol2./sum(reshape(vol5,N5,2),2)/4)-1) > 1e-6))
    %unique(C(tmp5(1:N5)))
    %unique(C(tmp5(N5+1:2*N5)))
    ItF = [true(N5,1); false(N5,1)]; ItL = not(ItF);
    tri2fac5t = tri2fac5'; fact = delfacc2(tri2fac5t);
    fact_ = tri2fac5t; fact_ = reshape(tri2fac5t(fact~=0),2,N5)';
    trinEA_ = reshape(fac2tri2(fact_),size(fact_)); %reshape when N5==1
    fact = reshape(fact(fact~=0),2,N5)'; qOm = fac2q(fact(:)); facO = fac2(fact(:),1:4);
    flp_ = trinEA_ == repmat(fI5,1,2);
    facO(flp_(:),:) = facO(flp_(:),[2 1 4 3]); qOm(flp_(:)) = not(qOm(flp_(:)));
    It_ = ndlockA(:,2) == fac2(fact(:,1),5); ItF_ = [It_ false(N5,1)]; ItL_ = [false(N5,1) It_];
    tmp = facO(ItF_,:); facO(ItF_,:) = facO(ItL_,:); facO(ItL_,:) = tmp;
    tmp = qOm(ItF_); qOm(ItF_) = qOm(ItL_); qOm(ItL_) = tmp;
    flp = facO(:,3) == [newtri5(ItF,3); newtri5(ItL,4)];
    facO(flp,:) = facO(flp,[2 1 4 3]); qOm(flp) = not(qOm(flp));
    facO = facO(:,[3 4 2 1]); nqOm = qOm; qOm = not(nqOm); %could be avoided
    
    %%facO_ =facO; qOm_ = qOm; %facO2 = facO_(:,[2 1 4 3]);
    %facO = [[newtri5(ItF,[4 3]); newtri5(ItL,[3 4])] [ndfreeA;ndfreeA]]; 
    %It_ = or(C==1,C==2); It = repmat(It_,2,1);
    %facO(It,[1 2 3 4]) = facO(It,[1 2 4 3]);
    %%It_ = any([C==1,C==3,C==5],2); ItF_ = [It_ false(N5,1)]; ItL_ = [false(N5,1) It_];
    %%tmp = facO(ItF_,:); facO(ItF_,:) = facO(ItL_,:); facO(ItL_,:) = tmp;
           %
    %%[tmp,Vo1a] = elem_inv(facO(:,[3 4 2]),xy);
    %%[tmp,Vo1b] = elem_inv(facO(:,[1 2 4]),xy); 
    %%[tmp,Vo1a] = elem_inv(facO(:,[4 3 1]),xy);
    %%[tmp,Vo1b] = elem_inv(facO(:,[2 1 3]),xy); 
    %%Vo1ab = Vo1a+Vo1b; Vo1ab = reshape(Vo1ab,N5,2);
    %%[tmp,area] = elem_inv([tri(I,[1 2 3]); tri(I,[1 3 4]); tri(I,[1 2 4]); tri(I,[2 3 4])],xy); area = sum(reshape(area,N5,4),2);
    %%[tmp,area] = elem_inv([[ndlockA(:,1); ndlockA(:,2)] facO(:,[3 4])],xy); area = reshape(area,N5,2);
    %%[sum(area,2)./sum(Vo1ab,2)*3/4 C] % [area./Vo1ab*3/4 C]
%
     %%figure(); ijk=24; lA = [facO(ijk,[3 2]) ndshr(ijk) facO(ijk,[1 4])]; lB = [facO(ijk+N5,[3 2]) ndshr(ijk) facO(ijk+N5,[1 4])]; plot3(xy(lA,1),xy(lA,2),xy(lA,3),'-b',xy(lB,1),xy(lB,2),xy(lB,3),'-r');
    %%qualO1 = elem_qual([facO(:,[3 2 4]); ...
    		      %%facO(:,[1 2 4])],xy,Nmetric,options); 
    %%qualO2 = elem_qual([facO(:,[4 1 3]); ...
    		      %%facO(:,[2 1 3])],xy,Nmetric,options);
    %%qualO1m = min(reshape(qualO1,N5*2,2),[],2);
    %%qualO2m = min(reshape(qualO2,N5*2,2),[],2);
    %%qOm = qualO1m > qualO2m; nqOm = not(qOm);
    %[qOm,nqOm] = surfqual(facO,[[3 2 4];[1 2 4]],[[4 1 3]; [2 1 3]],xy,Nmetric,options); 
    %ItF = [true(N5,1); false(N5,1)]; ItL = not(ItF);
    
    % none split to shared edge (next 25 % of cases)
    ItA = all(reshape(qOm ,N5,2),2);
    ItB = all(reshape(nqOm,N5,2),2);
    ItC = and(qOm(ItF),nqOm(ItL));
    ItD = and(nqOm(ItF),qOm(ItL));

    for i=1:4
    	switch i
    	case 1
    	pos = ItA; i1=4; j1=3;%all to facO(:,4)
    	case 2
    	pos = ItB; i1=3; j1=4;%all to facO(:,3)
    	case 3
    	pos = ItC; %first to facO(:,4)
    	otherwise
    	pos = ItD; %first to facO(:,3)
    	end;  
    	if not(any(pos))
    		continue;
    	end;
    N5_ = sum(pos); i2=i1-2;
    ItF_ = [pos; false(N5,1)]; ItL_ = [false(N5,1); pos];
    %note: facO(ItF,4) == facO(ItL,4) , facO(ItF,3) == facO(ItL,3) 
    TriT1 = [[facO(ItL_,[1 2]) facO(ItF_,1) ndshr(pos)]; ...
             [facO(ItF_,[2 1]) facO(ItL_,2) ndshr(pos)]]; %F1 to L2
    TriT2 = [[facO(ItL_,[1 2]) facO(ItF_,2) ndshr(pos)]; ...
             [facO(ItF_,[2 1]) facO(ItL_,1) ndshr(pos)]]; %F2 to L1
    if i==3
    TriT1 = [TriT2; [facO(ItF_,[2 1 4]) facO(ItL_,1)]; ...
                    [facO(ItL_,[1 2 3]) facO(ItF_,2)]; ...
                    [facO(ItF_,[4 2]) facO(ItL_,[1 3])]];
    elseif i==4
    TriT1 = [TriT1; [facO(ItL_,[1 2 4]) facO(ItF_,1)]; ...
                    [facO(ItF_,[2 1 3]) facO(ItL_,2)]; ...
                    [facO(ItL_,[2 4]) facO(ItF_,[1 3])]];    
    else
    TriT1 = [TriT1; [facO(ItL_,[1 2 i1]) facO(ItF_,1)]; ...
                    [facO(ItF_,[2 1 i1]) facO(ItL_,2)]; ...
                    [facO(ItF_,[4 3 i2]) facO(ItL_,i2)]];
    TriT2 = [TriT2; [facO(ItL_,[1 2 i1]) facO(ItF_,2)]; ...
                    [facO(ItF_,[2 1 i1]) facO(ItL_,1)]; ...
                    [facO(ItF_,[4 3 i2]) facO(ItL_,i2)]];
    end;
    %[tmp1,vol1] = elem_inv(TriT1,xy); vol1 = sum(reshape(vol1,sum(pos),5),2); %plot(tmp1(:))
    %[tmp2,vol2] = elem_inv(TriT2,xy); vol2 = sum(reshape(vol2,sum(pos),5),2);
    %[tmp,vol] = elem_inv(tri(I,:),xy);
    %[sum(vol(pos)) sum(vol1)*4/3 any(tmp1)]
    oldtri_ = repmat(fI5(pos),5,1);
    if i==1 || i == 2
     qual1 = elem_qual(TriT1,xy,Nmetric,options); 
     qual2 = elem_qual(TriT2,xy,Nmetric,options);
     qual1m = min(reshape(qual1,N5_,5),[],2);
     qual2m = min(reshape(qual2,N5_,5),[],2);
     option1 = qual1m > qual2m; option2 = not(option1); 
     option1 = repmat(option1,5,1);      
     option2 = repmat(option2,5,1); 
     q5 = [q5; qual1(option1); qual2(option2)];
     newtri5 = [newtri5; TriT1(option1,:); TriT2(option2,:)];
     %newtri5_ = [TriT1(option1,:); TriT2(option2,:)];
     oldtri = [oldtri; oldtri_(option1); oldtri_(option2)];     
    else
     q5 = [q5; elem_qual(TriT1,xy,Nmetric,options)];
     newtri5 = [newtri5; TriT1];
     %newtri5_ = TriT1;
     oldtri = [oldtri; oldtri_];
    end;
    end;%fori
    %[tmp,vol] = elem_inv(tri(I,:),xy);
    %[tmp,vol_] = elem_inv(newtri5,xy); [sum(vol) sum(vol_)]
    


function [newtri,q6] = mesh6long(tri_,ond,xy,Nmetric,options)
%nodes = [(1,3,6),(2,4,6),(3,2,5),(4,1,5)] in terms of edge numbers (magic of gen_books3D)
N6 = size(ond,1);
newtri = [[tri_(:,2) ond(:,[1 3 6])]; ...
	 [tri_(:,3) ond(:,[2 4 6])]; ...
	 [tri_(:,4) ond(:,[2 3 5])]; ...
	 [tri_(:,1) ond(:,[1 4 5])]];
%14 23 13 24 34 12: edges 1,3,5 opposes edges 2,4,6
ind = [[1 2 5 3];[1 2 3 6];[1 2 4 5];[1 2 6 4]];
I6 = [1 2 3 4 5 6; ...
      5 6 1 2 3 4; ...
      4 3 5 6 2 1];
newtri6 = zeros(4*N6,4); q6 = zeros(size(newtri,1),1); qmin = zeros(N6,1);
for i=1:3
	I6_ = I6(i,:);
	newtri6_ = reshape(ond(:,I6_(ind)),4*N6,4);
	q6_ = elem_qual(newtri6_,xy,Nmetric,options);
	qmin_ = min(reshape(q6_,N6,4),[],2);
	better = qmin < qmin_;
	qmin(better) = qmin_(better); 
	[tmp_,vol_] = elem_inv(newtri6_,xy);
	better6 = repmat(better,4,1);
	newtri6(better6,:) = newtri6_(better6,:);
	q6(better6) = q6_(better6);
end;%fori
q6 = [elem_qual(newtri,xy,Nmetric,options); q6];
newtri = [newtri; newtri6];
