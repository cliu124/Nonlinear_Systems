function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newfac,newedg] = fill_spheres(spheres,Nmetric,xy,badnds,badIDs,triQtb,options)
spheres = fix_spheres(spheres);
spheres = spheres(:,1:max(sum(squeeze(spheres(1,:,:))~=0)),:);
if options.prag_crs
	[newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_coarse(spheres,Nmetric,xy,badnds,badIDs,triQtb,options); return; end;
NN = sum(squeeze(spheres(1,:,:))~=0)'; 
badnds_start = badnds;
newtri = []; newtriN = []; qualityN = []; badbadnds = []; newIDs = [];
newedg = zeros(0,2); newedgN = []; newfac = zeros(0,3); newfacN = [];
if options.area == 0. %curved geometry
	[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds);
	%edg2tri_ = test_spheres(spheres,xy,badnds);
end;
while true
while true
[newtri,newfac_,spheres,badbadnds,newtriN,newIDs,newfacN_,qualityN,badnds,badIDs,NN,ndone] = rm_frnds(spheres,xy,NN,badnds,badIDs,badbadnds,Nmetric,triQtb,newtriN,newIDs,qualityN,newtri,options); %save for_debug3D.mat;
newfac = [newfac; newfac_];
newfacN = [newfacN; newfacN_];
%if ndone && size(spheres,3) ~= 0
%	spheres = fix_spheres(spheres);
%else 
if not(ndone) || size(spheres,3) == 0
	break;
end;
end; %while
if size(spheres,3) == 0 %nnz(spheres(1,:,:))
	break;
end;

%save for_debug3Db.mat;
%make an edge (and an element)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1) ~= edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg2tri = [tris(d) tris(ad)];
edg = edg(d,:);

Iflp = Ieven(d); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);
edgalt = [spheres(:,edg2tri(:,1)+(edg(:,1)-1)*size(spheres,2)); spheres(:,edg2tri(:,2)+(edg(:,1)-1)*size(spheres,2))];
edgalt(edgalt==repmat(edg(:,2)',6,1)) = 0;
edgalt(edgalt==repmat(edg(:,3)',6,1)) = 0;
edgalt = reshape(edgalt(edgalt~=0),2,size(edg,1))';
newtri_ = [edgalt edg(:,2:3)];
qual_ = elem_qual(newtri_,xy,Nmetric,options,[]);
if any(options.RMnd3D == [1 2])
	qual_(gen_badedgalt(edgalt,edg)) = 0; %we prevent splitting of spheres
end;


nN_ = [0; find([edg(1:end-1,1)~=edg(2:end,1); true])];
nN = diff(nN_); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN;
if options.qualRM
qual = zeros(numel(nN),nMax); qual(edg(:,1)+(C(I)-1)*numel(nN)) = qual_;
[maxqual,CI] = max(qual,[],2);
else
qual2_ = 1./elem_qual(edgalt,xy,Nmetric,options);
qual2_(qual_<=0) = 0; %no inversion
qual2 = zeros(numel(nN),nMax); qual2(edg(:,1)+(C(I)-1)*numel(nN)) = qual2_;
[maxqual,CI] = max(qual2,[],2);
end;
qualityN_ = zeros(size(CI)); Igood = maxqual>0;
if any(Igood)
qualityN_(Igood) = qual_(CI(Igood)+nN_([Igood; false]));
end;
if options.consRM
	Ibetter = qualityN_ > triQtb(badnds); %badnds(R2R)
else
	Ibetter = qualityN_ > options.minqual;
end;
if not(all(Ibetter))
badbadnds = [badbadnds; badnds(not(Ibetter))]; 
end;
if any(Ibetter)
RC = CI(Ibetter)+nN_([Ibetter; false]);
else
RC = [];
end;
newtri = [newtri; newtri_(RC,:)];
qualityN = [qualityN; qualityN_(Ibetter)];
newtriN = [newtriN; badnds(Ibetter)]; 
newIDs = [newIDs; badIDs(Ibetter)];
NN = NN(Ibetter);
badnds = badnds(Ibetter); 
badIDs = badIDs(Ibetter);
newedg = [newedg; newtri_(RC,[1 2])];
newedgN = [newedgN; badnds];
newfac = [newfac; newtri_(RC,[1 2 3]); newtri_(RC,[1 2 4])];
newfacN = [newfacN; repmat(badnds,2,1)];
spheres(:,edg2tri(RC,1)+(edg(RC,1)-1)*size(spheres,2)) = [edgalt(RC,:) edg(RC,2)]';
spheres(:,edg2tri(RC,2)+(edg(RC,1)-1)*size(spheres,2)) = [edgalt(RC,[2 1]) edg(RC,3)]';
spheres = spheres(:,:,Ibetter);
if size(spheres,3) == 0 
	break;
end;
if not(any(options.RMnd3D == [1 2]))
[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds);
end;
if options.debug == 2
edg2tri = test_spheres(spheres,xy,badnds);
end;
end; %while
Ikeep = true(max(badnds_start),1); 

Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs  = newIDs( Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
newfac = newfac(Ikeep(newfacN),:);
newfacN = newfacN(Ikeep(newfacN));
newedg = newedg(Ikeep(newedgN),:);
newedgN = newedgN(Ikeep(newedgN));

function [newtri,newfac,spheres,badbadnds,newtriN,newIDs,newfacN,qualityN,badnds,badIDs,NN,ndone] = rm_frnds(spheres,xy,NN,badnds,badIDs,badbadnds,Nmetric,triQtb,newtriN,newIDs,qualityN,newtri,options)
ndone = 0; newfac = zeros(0,3); newfacN = [];
R_ = reshape(repmat(1:size(spheres,3),3*size(spheres,2),1),size(spheres));
tris = reshape(repmat(1:size(spheres,2),3,size(spheres,3)),size(spheres));
Io = spheres~=0; tris = tris(Io);
nd2tri = [R_(Io) spheres(Io)];

[nd2tri,I] = sortrows(nd2tri); tris = tris(I);
I3 = and(and(nd2tri(2:end-1,1)==nd2tri(1:end-2,1), nd2tri(2:end-1,2)==nd2tri(1:end-2,2)), and(nd2tri(2:end-1,1)==nd2tri(3:end,1), nd2tri(2:end-1,2)==nd2tri(3:end,2)));
I3_ = [or(nd2tri(2,1)~=nd2tri(4,1),nd2tri(2,2)~=nd2tri(4,2)); ...
and(or(nd2tri(3:end-2,1)~=nd2tri(1:end-4,1), nd2tri(3:end-2,2)~=nd2tri(1:end-4,2)), or(nd2tri(3:end-2,1)~=nd2tri(5:end,1),nd2tri(3:end-2,2)~=nd2tri(5:end,2))); ...
or(nd2tri(end-1,1)~=nd2tri(end-3,1),nd2tri(end-1,2)~=nd2tri(end-3,2))];
d1 = and(I3,I3_);
if not(any(d1))
	return
end; 

%fix in case of tetrahedrons in spheres
I_ = NN==4;
if any(I_)
	[tmp,fC] = max(spheres(1,:,:)~=0,[],2); %no fix spheres
	first = reshape( repmat(spheres(1,squeeze(fC)+(0:size(spheres,3)-1)'*size(spheres,2)),3*size(spheres,2),1),size(spheres));
	%first = reshape( repmat(squeeze(spheres(1,1,:))',3*size(spheres,2),1),size(spheres));
	first = first(Io); first = first(I);
	I_ = and(I_(nd2tri(:,1)),nd2tri(:,2)~=first);
	d1(I_(2:end-1)) = false;
end;

d2 = [false; false; d1];
d3 = [d1; false; false];
d1 = [false; d1; false];
t1 = tris(d1)+(nd2tri(d1,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
t2 = tris(d2)+(nd2tri(d2,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
t3 = tris(d3)+(nd2tri(d3,1)-1)*size(spheres,2)*(size(spheres,3)~=1);
newtri1 = spheres(:,t1);
newtri2 = spheres(:,t2);
newtri3 = spheres(:,t3);
I1 = newtri1 == repmat(nd2tri(d1,2)',3,1); [C1,R] = find(I1);
I2 = newtri2 == repmat(nd2tri(d2,2)',3,1); [C2,R] = find(I2);
I3 = newtri3 == repmat(nd2tri(d3,2)',3,1); [C3,R] = find(I3);
thmp = [3 1 2]'; C1a = thmp(C1);
thmp = [2 3 1]'; C1 = thmp(C1); C2 = thmp(C2); C3 = thmp(C3);
R = (R-1)*3; newtri_ = [newtri1(C1+R) newtri2(C2+R) newtri3(C3+R) nd2tri(d1,2)];
Ieven = newtri1(C1a+R) ~= newtri2(C2+R);
newtri_(Ieven,1:2) = newtri_(Ieven,[2 1]); 
qualityN_ = elem_qual(newtri_,xy,Nmetric,options,[]);

Rs = nd2tri(d1,1); 
if numel(Rs) ~= 1
d = [Rs(1:end-1)~=Rs(2:end); true];
else
d = true;
end; %save for_debug3D.mat;
nN_ = [0; find(d)]; nN = diff(nN_); nMax = max(nN);
Rs2R = zeros(max(Rs),1); Rs2R(Rs(d)) = 1:numel(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
CI = reshape(C(I),numel(Rs),1); %reshape relevant for nMax==1
qual = inf(numel(nN),nMax); qual(Rs2R(Rs)+(CI-1)*numel(nN)) = qualityN_;
Rsd = Rs(d);

if options.consRM || nargin == 7
Ibetter = min(qual,[],2) > triQtb(badnds(Rsd));
else
Ibetter = min(qual,[],2) > options.minqual;
end;

if any(options.RMnd3D == [2 3]) %we will not kill spheres
IR = Ibetter(Rs2R(Rs)); nN = nN(Ibetter);
Rs = Rs(IR); d = d(IR);
t1 = t1(IR); t2 = t2(IR); t3 = t3(IR);
newtri_ = newtri_(IR,:);
qualityN_ = qualityN_(IR);
end;

%the following three lines are not dangerous, when you account for tetrahedrons in spheres
if any(d)
ndone = 1;
spheres(:,t1) = newtri_(:,1:3)';
spheres(:,t2) = 0;
spheres(:,t3) = 0;
newtriN_ = badnds(Rs);
newIDs_ = badIDs(Rs);
newtriN = [newtriN; newtriN_];
newIDs = [newIDs; newIDs_];
newtri = [newtri; newtri_];
qualityN = [qualityN; qualityN_];
I1 = or(NN(Rs)>6,NN(Rs)==4); I2 = NN(Rs)==6;
newfac = newtri_(I1,1:3); %incomplete
newfacN = newtriN_(I1); %incomplete
NN(Rs(d)) = NN(Rs(d))-2*nN;
end;
if any(not(Ibetter)) && not(any(options.RMnd3D == [2 3]))
Ibttr = true(size(spheres,3),1); Ibttr(Rsd(not(Ibetter))) = false;
badbadnds = [badbadnds; badnds(not(Ibttr))];
spheres = spheres(:,:,Ibttr);
badnds = badnds(Ibttr); 
badIDs = badIDs(Ibttr);
NN = NN(Ibttr);
end;
%finishes spheres
I = NN>3; spheres = spheres(:,:,I); NN = NN(I); badnds = badnds(I); badIDs = badIDs(I);
if size(spheres,3)~=0 && options.debug == 2
	edg2tri = test_spheres(spheres);
end;

function badedgalt = gen_badedgalt(edgalt,edg) %no semi intersecting spheres
nmax = max(max(edg(:,2:3)));
rmax = max(edg(:,1));
edg(:,1) = edg(:,1) + nmax;
nd2edg = inv_table(edg(:,1));
edgalt = sort(edgalt,2); 
inds = nd2edg(edg(:,1),:); I = inds ~= 0;
edg1 = zeros(size(inds)); edg1(I) = edg(inds(I),2);
edg2 = zeros(size(inds)); edg2(I) = edg(inds(I),3);
badedgalt = any(and(edg1==repmat(edgalt(:,1),1,size(nd2edg,2)),edg2==repmat(edgalt(:,2),1,size(nd2edg,2))),2);

%A1 = sparse(edg(:,2),edg(:,3),edg(:,1),nmax,nmax);
%A3 = sparse(edgalt(:,1),edgalt(:,2),1:size(edg,1),nmax,nmax);
%%edg and edgalt might be dubplicated across different sets!
%A2 = sparse(edgalt(:,1),edgalt(:,2),edg(:,1),nmax,nmax);
%I = and(A1==A2,A1~=0);
%badedgalt = full(A3(I));
%

function ospheres = fix_spheres(spheres)
spheres = permute(spheres,[2 1 3]);
[Cg,Rg] = find(spheres~=0);
[C ,R ] = find(spheres==0);
move = zeros(size(spheres)); move(C+(R-1)*size(spheres,1)) = -1; move = cumsum(move);
Cg = Cg + reshape(move(Cg+(Rg-1)*size(spheres,1)),size(Cg));
ospheres = zeros(size(spheres));
ospheres(Cg+(Rg-1)*size(ospheres,1)) = spheres(spheres~=0);
ospheres = permute(ospheres,[2 1 3]);

function [spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,NN,badnds,badIDs,badbadnds)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
edg = sort(edg,2);
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); 
d = [false; false; and(and(edg(3:end-1,1)==edg(4:end,1),edg(3:end-1,1)==edg(1:end-3,1)), ...
and(and(edg(3:end-1,2)==edg(4:end,2),edg(3:end-1,2)==edg(1:end-3,2)), ...
    and(edg(3:end-1,3)==edg(4:end,3),edg(3:end-1,3)==edg(1:end-3,3)))); false];
if any(d)
	Rbad = edg(d,1);
	badbadnds = [badbadnds; badnds(Rbad)]; 
	Ikeep = true(size(spheres,3),1); Ikeep(Rbad) = false;
	spheres = spheres(:,:,Ikeep); badnds = badnds(Ikeep); badIDs = badIDs(Ikeep);
	if numel(NN) ~= 0
	NN = NN(Ikeep);  
	end;
end;

function edg2tri = test_spheres(spheres,xy,badnds)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
fac = sortrows([R_(I) sort([spheres(1,I)' spheres(2,I)' spheres(3,I)'],2)]);
d = all(fac(1:end-1,:) == fac(2:end,:),2);
if any(d)
	save for_debug3D.mat; 
	error('double face in spheres');
end;
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1)~=edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg = edg(d,:);
edg2tri = [tris(d) tris(ad)];

function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_coarse(spheres,Nmetric,xy,badnds,badIDs,triQtb,options)
badbadnds = zeros(0,1); 
if options.area == 0. %curved geometry
	[spheres,NN,badnds,badIDs,badbadnds] = kill_split_spheres(spheres,[],badnds,badIDs,badbadnds);
	%edg2tri_ = test_spheres(spheres,xy,badnds);
end;
newtri = zeros(0,4); newtriN = []; qualityN = []; badnds_start = badnds; newIDs = [];

nbad = size(spheres,3); nbad_ = size(spheres,2); nbadE = 3*nbad_;
[Cc,R] = find(reshape(spheres,[nbadE,size(spheres,3)])~=0); 
edg = [badnds(R) reshape(spheres(Cc+(R-1)*nbadE),size(R))];
[edg_,I] = sortrows(edg); [tmp,Iback] = sort(I);
Igood = [true; or(edg_(1:end-1,1) ~= edg_(2:end,1),edg_(1:end-1,2) ~= edg_(2:end,2))];
edg = edg(Igood(Iback),:); R = R(Igood(Iback)); Cc = Cc(Igood(Iback));
edgL = inf(  nbad,nbadE);
edgS = zeros(nbad,nbadE);
edgL(R+(Cc-1)*nbad) = elem_qual(edg,xy,Nmetric,options);
edgS(R+(Cc-1)*nbad) = edg(:,2);
[edgL_,Iedg] = sort(edgL,2);
edgS = reshape(edgS(reshape(repmat((1:nbad)',1,nbadE),nbad*nbadE,1)+(Iedg(:)-1)*nbad),size(edgS));
%Iedg(isinf(edgL)) = 0; 
theC = 1;
NN = sum(edgS~=0,2);
while true
[Cc,R] = find(squeeze(spheres(1,:,:))~=0); nbad = size(spheres,3);
if nbad==1
	R_ = R; R = reshape(Cc,numel(Cc),1); Cc = reshape(R_,size(R));
end;
%C_ = Iedg(R+(theC-1)*nbad);
thshrt = reshape(edgS(R+(theC-1)*nbad),numel(R),1);
newtri_ = [thshrt spheres(1,Cc+(R-1)*nbad_)' spheres(2,Cc+(R-1)*nbad_)' spheres(3,Cc+(R-1)*nbad_)'];
[newtri__,tmp] = sort(newtri_,2);
Igood = and(and(newtri__(:,1)~=newtri__(:,2),newtri__(:,2)~=newtri__(:,3)),newtri__(:,3)~=newtri__(:,4));
newtri_ = newtri_(Igood,:); Cc = Cc(Igood); R = R(Igood);
newqual = elem_qual(newtri_,xy,Nmetric,options,1);
if options.consRM
	Igood = newqual > triQtb(badnds(R));
else
	Igood = newqual > options.minqual;
end;
allgood = true(nbad,nbad_); allgood(R+(Cc-1)*nbad) = Igood;
Igood = all(allgood,2); Ibad = not(Igood);
newtri = [newtri; newtri_(Igood(R),:)];
qualityN = [qualityN; newqual(Igood(R))];
newtriN  = [newtriN ; badnds(R(Igood(R)))];
newIDs   = [newIDs  ; badIDs(R(Igood(R)))];
%update circles
spheres = spheres(:,:,Ibad); badnds = badnds(Ibad); NN = NN(Ibad); edgS = edgS(Ibad,:); badIDs = badIDs(Ibad);
if size(spheres,3) == 0
	break;
end;
theC = theC+1;
Ibad = NN < theC; Igood = not(Ibad);
badbadnds = [badbadnds; badnds(Ibad)];
spheres = spheres(:,:,Igood); badnds = badnds(Igood); NN = NN(Igood); edgS = edgS(Igood,:); badIDs = badIDs(Igood);
if size(spheres,3) == 0
	break;
end; %[size(spheres,3) size(edgS,1) theC]
end;
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
%next four lines are irrelevant outside the context of region IDs
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));