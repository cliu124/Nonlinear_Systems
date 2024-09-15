function [newtri,qualityN,badnds,newtriN,badbadnds,newfac,newedg] = fill_spheres_new(spheres,Nmetric,xy,badnds,triQtb,options)
badnds_start = badnds;
spheres = fix_spheres(spheres);
spheres = spheres(:,1:max(sum(squeeze(spheres(1,:,:))~=0)),:);
NN = sum(squeeze(spheres(1,:,:))~=0)'; %R2R = (1:numel(badnds))';
newtri = []; newtriN = []; qualityN = []; badbadnds = []; %R2R = 1:numel(badnds);
newedg = zeros(0,2); newedgN = []; newfac = zeros(0,3); newfacN = [];

%make an edge (and an element)
I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1)~=edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg2tri = [tris(d) tris(ad)];
edg = edg(d,:);

Iflp = Ieven(d); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);
edgalt = [spheres(:,edg2tri(:,1)+(edg(:,1)-1)*size(spheres,2)); spheres(:,edg2tri(:,2)+(edg(:,1)-1)*size(spheres,2))];
edgalt(edgalt==repmat(edg(:,2)',6,1)) = 0;
edgalt(edgalt==repmat(edg(:,3)',6,1)) = 0;
edgalt = reshape(edgalt(edgalt~=0),2,size(edg,1))';
swaptri = [edgalt edg];
swapqual = elem_qual(swaptri(:,[1 2 4 5]),xy,Nmetric,options,[]);

while true
while true
[newtri,newfac_,spheres,badbadnds,newtriN,newfacN_,qualityN,badnds,NN,edg2tri,swaptri,swapqual,ndone] = rm_frnds(spheres,xy,NN,badnds,badbadnds,Nmetric,triQtb,newtriN,qualityN,newtri,edg2tri,swaptri,swapqual,options); 
%save for_debug3D.mat; error('just stop3D');
if options.debug
test_swap(swaptri,spheres,swapqual,edg2tri,xy,Nmetric,options);
end;
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

if any(options.RMnd3D == [1 2])
	swapqual(gen_badedgalt(swaptri(:,[1 2]),swaptri(:,[3 4 5]))) = -1; %we prevent splitting of spheres
end;

swapqual(isnan(swapqual)) = elem_qual(swaptri(isnan(swapqual),[1 2 4 5]),xy,Nmetric,options,[]);
nN_ = [0; find([swaptri(1:end-1,3)~=swaptri(2:end,3); true])];
nN = diff(nN_); nMax = max(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN;
if options.qualRM
qual = zeros(numel(nN),nMax); qual(swaptri(:,3)+(C(I)-1)*numel(nN)) = swapqual;
[maxqual,CI] = max(qual,[],2);
else
qual2_ = 1./elem_qual(swaptri(:,[1 2]),xy,Nmetric,options);
qual2_(swapqual<=0) = 0; %no inversion
qual2 = zeros(numel(nN),nMax); qual2(swaptri(:,3)+(C(I)-1)*numel(nN)) = qual2_;
[maxqual,CI] = max(qual2,[],2);
end;
qualityN_ = zeros(size(CI)); Igood = maxqual>0;
if any(Igood)
qualityN_(Igood) = swapqual(CI(Igood)+nN_([Igood; false]));
end;
if options.consRM
	Ibetter = qualityN_ > triQtb(badnds); %badnds(R2R)
else
	Ibetter = qualityN_ > options.minqual;
end;
if not(all(Ibetter))
badbadnds = [badbadnds; badnds(not(Ibetter))]; %badnds(R2R(not(Ibetter)))
end;
if any(Ibetter)
RC = CI(Ibetter)+nN_([Ibetter; false]);
else
RC = [];
end;
newtri = [newtri; swaptri(RC,[1 2 4 5])];
qualityN = [qualityN; qualityN_(Ibetter)];
newtriN = [newtriN; badnds(Ibetter)]; %badnds(R2R(Ibetter))

NN = NN(Ibetter);
badnds = badnds(Ibetter); %R2R = R2R(Ibetter);
newedg = [newedg; swaptri(RC,[1 2])];
newedgN = [newedgN; badnds];
newfac = [newfac; swaptri(RC,[1 2 4]); swaptri(RC,[1 2 5])];
newfacN = [newfacN; repmat(badnds,2,1)];
spheres(:,edg2tri(RC,1)+(swaptri(RC,3)-1)*size(spheres,2)) = [swaptri(RC,[1 2]) swaptri(RC,4)]';
spheres(:,edg2tri(RC,2)+(swaptri(RC,3)-1)*size(spheres,2)) = [swaptri(RC,[2 1]) swaptri(RC,5)]';
[swaptri,swapqual,edg2tri] = fix_edg(swaptri,swapqual,edg2tri,RC);


spheres = spheres(:,:,Ibetter);
I = Ibetter(swaptri(:,3));
swaptri = swaptri(I,:); swapqual = swapqual(I); edg2tri = edg2tri(I,:);
old2new = zeros(size(Ibetter)); old2new(Ibetter) = 1:numel(NN); swaptri(:,3) = old2new(swaptri(:,3));
%old2new = zeros(size(I)); old2new(I) = 1:size(swaptri,1);
%nd2tri = nd2tri(Ibetter,:);
%nd2tri(nd2tri~=0) = old2new(nd2tri(nd2tri~=0));
%nd2tri = fix_circles(nd2tri);
if size(spheres,3) == 0 %nnz(spheres(1,:,:))
	break;
end;
if not(any(options.RMnd3D == [1 2]))
[spheres,NN,badnds,badbadnds] = kill_split_spheres(spheres,NN,badnds,badbadnds);
end;

if options.debug == 2
edg2tri_ = test_spheres(spheres);
end;
end; %while
Ikeep = true(size(xy,1),1); 

Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
newfac = newfac(Ikeep(newfacN),:);
newfacN = newfacN(Ikeep(newfacN));
newedg = newedg(Ikeep(newedgN),:);
newedgN = newedgN(Ikeep(newedgN));



function [newtri,newfac,spheres,badbadnds,newtriN,newfacN,qualityN,badnds,NN,edg2tri,swaptri,swapqual,ndone] = rm_frnds(spheres,xy,NN,badnds,badbadnds,Nmetric,triQtb,newtriN,qualityN,newtri,edg2tri,swaptri,swapqual,options)
ndone = 0; newfac = zeros(0,3); newfacN = [];
if options.debug
	test_swap(swaptri,spheres,swapqual,edg2tri,xy,Nmetric,options);
end;
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

Rs = nd2tri(d1,1); trisd1 = tris(d1);
if numel(Rs) ~= 1
d = [Rs(1:end-1)~=Rs(2:end); true];
else
d = true;
end; 
nN_ = [0; find(d)]; nN = diff(nN_); nMax = max(nN);
Rs2R = zeros(max(Rs),1); Rs2R(Rs(d)) = 1:numel(nN);
C = repmat((1:nMax)',1,numel(nN)); CnN = repmat(nN',nMax,1); I = C <= CnN; 
CI = reshape(C(I),numel(Rs),1); %reshape relevant for nMax==1
qual = inf(numel(nN),nMax); qual(Rs2R(Rs)+(CI-1)*numel(nN)) = qualityN_;
Rsd = Rs(d); 

if options.consRM || nargin == 7
Ibetter = min(qual,[],2) > triQtb(badnds(Rsd)); %badnds(R2R(Rsd))
else
Ibetter = min(qual,[],2) > options.minqual;
end;

if any(options.RMnd3D == [2 3]) %we will not kill spheres
IR = Ibetter(Rs2R(Rs)); nN = nN(Ibetter);
Rs = Rs(IR); d = d(IR); trisd1 = trisd1(IR); 
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
[swaptri,swapqual,edg2tri] = fix_edg2(swaptri,swapqual,Rs,newtri_,edg2tri,trisd1);
newtriN = [newtriN; newtriN_];
newtri = [newtri; newtri_];
qualityN = [qualityN; qualityN_];

I1 = or(NN(Rs)>6,NN(Rs)==4); I2 = NN(Rs)==6;
newfac = newtri_(I1,1:3); %incomplete in case of double tetrahedrons
newfacN = newtriN_(I1); %incomplete in case of double tetrahedrons
NN(Rs(d)) = NN(Rs(d))-2*nN;
end;
if any(not(Ibetter)) && not(any(options.RMnd3D == [2 3]))
Ibttr = true(size(spheres,3),1); Ibttr(Rsd(not(Ibetter))) = false;
badbadnds = [badbadnds; badnds(not(Ibttr))]; %badnds(R2R(not(Ibttr)))
spheres = spheres(:,:,Ibttr);
badnds = badnds(Ibttr); %R2R = R2R(Ibttr);
NN = NN(Ibttr);
I = Ibttr(swaptri(:,3));
swapqual = swapqual(I); swaptri = swaptri(I,:); edg2tri = edg2tri(I,:);
oR2nR = zeros(size(Ibttr)); oR2nR(Ibttr) = 1:numel(NN);
swaptri(:,3) = oR2nR(swaptri(:,3));
%old2new = zeros(size(I)); old2new(I) = 1:size(swaptri,1);
%nd2tri = nd2tri(Ibttr,:);
%nd2tri(nd2tri~=0) = old2new(nd2tri(nd2tri~=0));
%nd2tri = fix_circles(nd2tri);

end;
%finishes spheres
Irmn = NN>3; spheres = spheres(:,:,Irmn); NN = NN(Irmn); badnds = badnds(Irmn); %R2R = R2R(I);
I = Irmn(swaptri(:,3));
swapqual = swapqual(I); swaptri = swaptri(I,:); edg2tri = edg2tri(I,:);
%old2new = zeros(size(I)); old2new(I) = 1:size(swaptri,1);
%nd2tri = nd2tri(Ibttr,:);
%nd2tri(nd2tri~=0) = old2new(nd2tri(nd2tri~=0));
%nd2tri = fix_circles(nd2tri);
oR2nR = zeros(size(Irmn)); oR2nR(Irmn) = 1:numel(NN);
swaptri(:,3) = oR2nR(swaptri(:,3));
if size(spheres,3)~=0 && options.debug == 2
	edg2tri = test_spheres(spheres);
end;

function [swaptri,swapqual,edg2tri,nd2tri] = fix_edg2(swaptri,swapqual,Rs,newtri_,edg2tri,t1)
%save for_debug3D.mat; error('just stop3D');
nN_ = [0; find([swaptri(1:end-1,3)~=swaptri(2:end,3); true])];
IR = false(swaptri(end,3),1); IR(Rs) = true; IR = IR(swaptri(:,3)); 
swaptri_ = swaptri(IR,:);
nd2tri = inv_table(swaptri(:,3));
glb2lcl = zeros(size(swaptri,1),1); glb2lcl(IR) = 1:size(swaptri_,1);
nd2tri_ = nd2tri; nd2tri_(nd2tri_~=0) = glb2lcl(nd2tri_(nd2tri_~=0));
nd2tri_ = nd2tri_(:,1:max(sum(nd2tri_~=0,2)));
%nd2tri_ = inv_table(swaptri_(:,3));
inds = nd2tri_(Rs,:)'; I = inds ~= 0;
edg1 = zeros(size(inds)); edg1(I) = swaptri_(inds(I),1);
edg2 = zeros(size(inds)); edg2(I) = swaptri_(inds(I),2);
edg3 = zeros(size(inds)); edg3(I) = swaptri_(inds(I),4);
edg4 = zeros(size(inds)); edg4(I) = swaptri_(inds(I),5);
tri1 = repmat(newtri_(:,1)',size(I,1),1); 
tri2 = repmat(newtri_(:,2)',size(I,1),1);
tri3 = repmat(newtri_(:,3)',size(I,1),1);
tri4 = repmat(newtri_(:,4)',size(I,1),1);
%update swaptri, swapqual and edg2tri
I1a = and(and(edg3 == tri2,edg4 == tri3),edg2==tri4);
I1b = and(and(edg4 == tri2,edg3 == tri3),edg1==tri4);
I2a = and(and(edg3 == tri3,edg4 == tri1),edg2==tri4);
I2b = and(and(edg4 == tri3,edg3 == tri1),edg1==tri4);
I3a = and(and(edg3 == tri1,edg4 == tri2),edg2==tri4);
I3b = and(and(edg4 == tri1,edg3 == tri2),edg1==tri4);
%I1 = or(I1a,I1b); I2 = or(I2a,I2b); I3 = or(I3a,I3b);
%Ia = or(or(I1a,I2a),I3a); Ib = or(or(I1b,I2b),I3b);
%edg_ = [edg1(I1a) edg2(I1a) Rs(any(I1a)) edg3(I1a) edg4(I1a)];
I = any(I1a);
if any(I)
[C,R] = find(I1a); C = C+nN_(Rs(I));
swaptri(C,2) = tri1(I1a); swapqual(C) = NaN;
edg2tri(C,2) = t1(I);
end;
I = any(I1b);
if any(I)
[C,R] = find(I1b); C = C+nN_(Rs(I));
swaptri(C,1) = tri1(I1b); swapqual(C) = NaN;
edg2tri(C,1) = t1(I);
end;
I = any(I2a);
if any(I)
[C,R] = find(I2a); C = C+nN_(Rs(I));
swaptri(C,2) = tri2(I2a); swapqual(C) = NaN;
edg2tri(C,2) = t1(I);
end;
I = any(I2b);
if any(I)
[C,R] = find(I2b); C = C+nN_(Rs(I));
swaptri(C,1) = tri2(I2b); swapqual(C) = NaN;
edg2tri(C,1) = t1(I);
end;
I = any(I3a);
if any(I)
[C,R] = find(I3a); C = C+nN_(Rs(I));
swaptri(C,2) = tri3(I3a); swapqual(C) = NaN;
edg2tri(C,2) = t1(I);
end;
I = any(I3b);
if any(I)
[C,R] = find(I3b); C = C+nN_(Rs(I));
swaptri(C,1) = tri3(I3b); swapqual(C) = NaN;
edg2tri(C,1) = t1(I);
end;

%take out edges related to tri4 in swapqual
I = or(edg3 == tri4, edg4 == tri4);
[C,R] = find(I); C = C+nN_(Rs(R));
I = true(size(swaptri,1),1); I(C) = false;
swaptri = swaptri(I,:); swapqual = swapqual(I); edg2tri = edg2tri(I,:);
%fix nd2tri
old2new = zeros(size(I)); old2new(I) = 1:size(swaptri,1);
nd2tri(nd2tri~=0) = old2new(nd2tri(nd2tri~=0));
nd2tri = fix_circles(nd2tri);
%I = nd2tri == inv_table(swaptri(:,3)); sum(I(:)) == numel(nd2tri)


function [swaptri,swapqual,edg2tri] = fix_edg(swaptri,swapqual,edg2tri,RC)
%save for_debug3D.mat; error('just stop3D');
swapqual(RC) = -1;
I = swaptri(RC,1) < swaptri(RC,2); nI = not(I);
swaptri(RC(I) ,[1 2 4 5]) = swaptri(RC(I) ,[5 4 1 2]);
swaptri(RC(nI),[1 2 4 5]) = swaptri(RC(nI),[4 5 2 1]);
%swaptri(RC,1) was contained in spheres(:,edg2tri(RC,1)+(swaptri(RC,3)-1)*size(spheres,2))'
edg2tri(RC(I),:) = edg2tri(RC(I),[2 1]);
%it does not end here!
nd2tri = inv_table(swaptri(:,3));
%nd2tri_ = inv_table(swaptri_(:,3));
Rs = swaptri(RC,3); 
nN_ = [0; find([swaptri(1:end-1,3)~=swaptri(2:end,3); true])];
inds = nd2tri(Rs,:)'; I = inds ~= 0;
edg1 = zeros(size(inds)); edg1(I) = swaptri(inds(I),4);
edg2 = zeros(size(inds)); edg2(I) = swaptri(inds(I),5);
edg3 = zeros(size(inds)); edg3(I) = swaptri(inds(I),1);
edg4 = zeros(size(inds)); edg4(I) = swaptri(inds(I),2);
tri1 = repmat(swaptri(RC,4)',size(I,1),1); 
tri2 = repmat(swaptri(RC,5)',size(I,1),1);
tri3 = repmat(swaptri(RC,1)',size(I,1),1);
tri4 = repmat(swaptri(RC,2)',size(I,1),1);
%update swaptri, swapqual and edg2tri
I1a = and(and(edg1 == tri3,edg3 == tri4),edg2==tri1); %sum(I1a(:))
I1b = and(and(edg3 == tri3,edg1 == tri4),edg2==tri2); %sum(I1b(:))
I2a = and(and(edg1 == tri3,edg4 == tri4),edg2==tri2); %sum(I2a(:))
I2b = and(and(edg4 == tri3,edg1 == tri4),edg2==tri1); %sum(I2b(:))
I3a = and(and(edg2 == tri3,edg3 == tri4),edg1==tri2); %sum(I3a(:))
I3b = and(and(edg3 == tri3,edg2 == tri4),edg1==tri1); %sum(I3b(:))
I4a = and(and(edg2 == tri3,edg4 == tri4),edg1==tri1); %sum(I4a(:))
I4b = and(and(edg4 == tri3,edg2 == tri4),edg1==tri2); %sum(I4b(:))
%I1a = and(and(edg1 == tri3,edg3 == tri4),edg2==tri2); sum(I1a(:))
%I1b = and(and(edg3 == tri3,edg1 == tri4),edg2==tri1); sum(I1b(:))
%I2a = and(and(edg1 == tri3,edg4 == tri4),edg2==tri1); sum(I2a(:))
%I2b = and(and(edg4 == tri3,edg1 == tri4),edg2==tri2); sum(I2b(:))
%I3a = and(and(edg2 == tri3,edg3 == tri4),edg1==tri1); sum(I3a(:))
%I3b = and(and(edg3 == tri3,edg2 == tri4),edg1==tri2); sum(I3b(:))
%I4a = and(and(edg2 == tri3,edg4 == tri4),edg1==tri2); sum(I4a(:))
%I4b = and(and(edg4 == tri3,edg2 == tri4),edg1==tri1); sum(I4b(:))
I1 = or(I1a,I1b); I2 = or(I2a,I2b); I3 = or(I3a,I3b); I4 = or(I4a,I4b);
%swaptri(RC(1:10),:)
I = any(I1a);
if any(I)
[C,R] = find(I1a); C = C+nN_(Rs(I));
swaptri(C,1) = tri2(I1a); swapqual(C) = NaN;
edg2tri(C,1) = edg2tri(RC(I),1);
end;
I = any(I1b);
if any(I)
[C,R] = find(I1b); C = C+nN_(Rs(I));
swaptri(C,1) = tri1(I1b); swapqual(C) = NaN;
edg2tri(C,1) = edg2tri(RC(I),2);
end;
I = any(I2a);
if any(I)
[C,R] = find(I2a); C = C+nN_(Rs(I));
swaptri(C,2) = tri1(I2a); swapqual(C) = NaN;
%edg2tri(C,2) = edg2tri(RC(I),2); already ==
end;
I = any(I2b);
if any(I)
[C,R] = find(I2b); C = C+nN_(Rs(I));
swaptri(C,2) = tri2(I2b); swapqual(C) = NaN;
%edg2tri(C,2) == edg2tri(RC(I),1); already ==
end;
I = any(I3a);
if any(I)
[C,R] = find(I3a); C = C+nN_(Rs(I));
swaptri(C,1) = tri1(I3a); swapqual(C) = NaN;
%edg2tri(C,1) = edg2tri(RC(I),2); %already ==
end;
I = any(I3b);
if any(I)
[C,R] = find(I3b); C = C+nN_(Rs(I));
swaptri(C,1) = tri2(I3b); swapqual(C) = NaN;
%edg2tri(C,1) = edg2tri(RC(I),1); %already ==
end;
I = any(I4a);
if any(I)
[C,R] = find(I4a); C = C+nN_(Rs(I));
swaptri(C,2) = tri2(I4a); swapqual(C) = NaN;
edg2tri(C,2) = edg2tri(RC(I),1);
end;
I = any(I4b);
if any(I)
[C,R] = find(I4b); C = C+nN_(Rs(I));
swaptri(C,2) = tri1(I4b); swapqual(C) = NaN;
edg2tri(C,2) = edg2tri(RC(I),2);
end;
%edg2tri(RC,:) = edg2tri(RC,[2 1]);



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

function [spheres,NN,badnds,badbadnds] = kill_split_spheres(spheres,NN,badnds,badbadnds)
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
	Rbad = edg(d,1); %save for_debug3D.mat;
	%if numel(badnds) == 1
	badbadnds = [badbadnds; badnds(Rbad)]; %R2R(Rbad)
	%else
	%badbadnds = [badbadnds; badnds(Rbad)]; %R2R(Rbad)
	%end;
	Ikeep = true(size(spheres,3),1); Ikeep(Rbad) = false;
	spheres = spheres(:,:,Ikeep); NN = NN(Ikeep); badnds = badnds(Ikeep); %R2R = R2R(Ikeep);
	qual_ = qual_(Ikeep(newtri_(:,3)));
	newtri_ = newtri_(Ikeep(newtri_(:,3)),:);
	old2new = zeros(size(Ikeep)); old2new(Ikeep) = 1:nnz(Ikeep); newtri_(:,3) = old2new(newtri_(:,3));
end;

function edg2tri = test_spheres(spheres)
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

function [] = test_swap(swaptri,spheres,swapqual,edg2tri_,xy,Nmetric,options)
if not(all(swaptri(:,4) < swaptri(:,5)))
	error('largest node should be last in edge list');
end;
[tmp,I_] = sortrows(swaptri(:,3:5));
swaptri = swaptri(I_,:); swapqual = swapqual(I_); edg2tri_ = edg2tri_(I_,:);

I = squeeze(spheres(1,:,:))~=0;
R_ = repmat(1:size(spheres,3),size(spheres,2),1);
tris = reshape(repmat(1:size(spheres,2),1,size(spheres,3)),[size(spheres,2) size(spheres,3)]); tris = repmat(tris(I),3,1);
edg = [[spheres(1,I)'; spheres(2,I)'; spheres(3,I)'] [spheres(2,I)'; spheres(3,I)'; spheres(1,I)']];
[edg,Ieven] = sort(edg,2); Ieven = Ieven(:,1) == 1;
edg = [repmat(R_(I),3,1) edg];
[edg,I] = sortrows(edg); tris = tris(I); Ieven = Ieven(I);
d = [true; or(edg(1:end-1,1)~=edg(2:end,1),or(edg(1:end-1,2) ~= edg(2:end,2),edg(1:end-1,3) ~= edg(2:end,3)))]; ad = not(d);
edg2tri = [tris(d) tris(ad)];
edg = edg(d,:);

Iflp = Ieven(d); edg2tri(Iflp,:) = edg2tri(Iflp,[2 1]);
edgalt = [spheres(:,edg2tri(:,1)+(edg(:,1)-1)*size(spheres,2)); spheres(:,edg2tri(:,2)+(edg(:,1)-1)*size(spheres,2))];
edgalt(edgalt==repmat(edg(:,2)',6,1)) = 0;
edgalt(edgalt==repmat(edg(:,3)',6,1)) = 0;
edgalt = reshape(edgalt(edgalt~=0),2,size(edg,1))';
if any(any(edg ~= swaptri(:,3:5),2))
	error(sprintf('edg is flawed(%0.0f)',find(any(edg ~= swaptri(:,3:5),2),1)));
end;
if any(any(edgalt ~= swaptri(:,1:2),2))
	save for_debug3D.mat;
	error(sprintf('%0.0f edgalt is flawed(%0.0f)',sum(any(edgalt ~= swaptri(:,1:2),2)),find(any(edgalt ~= swaptri(:,1:2),2),1)));
end;
swapqual_ = elem_qual([edgalt edg(:,[2 3])],xy,Nmetric,options,[]);
swapqual(isnan(swapqual)) = elem_qual(swaptri(isnan(swapqual),[1 2 4 5]),xy,Nmetric,options,[]);
if sum(abs(swapqual_ -swapqual)>1e-6)
	error(sprintf('swapqual is flawed (%0.0f)',find(abs(swapqual_ -swapqual)>1e-6,1)));
end;
if any(any(edg2tri ~= edg2tri_,2))
	error(sprintf('%0.0f edg2tri is flawed(%0.0f)',sum(any(edg2tri ~= edg2tri_,2)),find(any(edg2tri ~= edg2tri_,2),1)));
end;



function ocircles = fix_circles(circles)
[C,R] = find(circles'~=0);
[Cg,Rg] = find((repmat(sum(circles~=0,2),1,size(circles,2)) >= repmat(1:size(circles,2),size(circles,1),1))');
ocircles = zeros(size(circles));
ocircles(Rg+(Cg-1)*size(circles,1)) = circles(R+(C-1)*size(circles,1));