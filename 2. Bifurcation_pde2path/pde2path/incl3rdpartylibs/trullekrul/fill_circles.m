function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs,newedg,newedgN,newfac,newfacN] = fill_circles(circles,Nmetric,xy,badnds,badIDs,triQtb,options,edg_)
NN = sum(circles~=0,2);
if size(edg_,2) == 3
circles = fix_circles(circles);
end;
if options.prag_crs && size(edg_,2) ~= 2
[newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_coarse(circles,Nmetric,xy,badnds,badIDs,triQtb,NN,options,edg_); return; end;
%if options.rmnonconvec~=2
%[circles,badnds,edg_,badbadnds] = del_dblNcnvx(circles,xy,badnds,edg_,options);
%else
badbadnds = zeros(0,1);
%end;
newtri = zeros(0,3); newtriN = []; qualityN = []; badnds_start = badnds;
newedg = zeros(0,2); newedgN = []; newIDs = [];
%[newtri,circles,badnds,badbadnds,newtriN,qualityN,R2R,NN,edg_] = rm_trinds(circles,badnds,badbadnds,xy,Nmetric,triQtb,options,edg_);
%if size(circles,1) > 0
%NN = sum(circles~=0,2);
qual_ = zeros(size(circles,2),size(circles,1)); qual_(circles'~=0) = NaN;
if size(edg_,2) == 2
newtri = zeros(0,4); 
newfac = zeros(0,3);
newfacN = [];
qual3Da = qual_;
qual3Db = qual_;
elseif not(options.qualRM) 
	qual2 = qual_;
end;
allconvex = true;
while true 
    % FIND ALLOWED NEW ELEMENTS
    [Cc,R] = find(circles'); nbad = size(circles,1);
    thmp1 = [0:size(circles,2)-1]'; thmp2 = [2:size(circles,2)+1]';
    Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
    Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
    I1 = circles(R+(Cr-1)*nbad); I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
    Ipdt = isnan(qual_(Cc+(R-1)*size(qual_,1))); I123 = [I1(:) I2(:) I3(:)];
    if size(xy,2) == 2
    qual_(isnan(qual_(:))) = elem_qual(I123(Ipdt,:),xy,Nmetric,options,1);
    elseif size(edg_,2) == 2 %swapping
    %qual3Da(isnan(qual_(:))) = elem_qual([edg_(R(Ipdt),1) I123(Ipdt,[3 2 1])], xy,Nmetric,options,1);
    %qual3Db(isnan(qual_(:))) = elem_qual([edg_(R(Ipdt),2) I123(Ipdt,[2 3 1])], xy,Nmetric,options,1);
    qual3Da(isnan(qual_(:))) = elem_qual([I123(Ipdt,[1 2 3]) edg_(R(Ipdt),1)], xy,Nmetric,options,1);
    qual3Db(isnan(qual_(:))) = elem_qual([I123(Ipdt,[3 2 1]) edg_(R(Ipdt),2)], xy,Nmetric,options,1);
    qual3Dmin = min(qual3Da,qual3Db);
    qual_(isnan(qual_(:))) = qual3Dmin(isnan(qual_(:)));
    else %3D
    qual_(isnan(qual_(:))) = elem_qual(I123(Ipdt,:),xy,Nmetric,options,edg_(R(Ipdt),:));
    end;
    
    if options.qualRM || size(edg_,2) == 2
	qual_fix = qual_;
    else 
    	qual2(isnan(qual2(:))) = 1./elem_qual(I123(Ipdt,[1 3]),xy,Nmetric,options);
    	qual_fix = qual2;
    	qual_fix(qual_ < 0) = 0;
    end;
    if not(allconvex) %only triangulate next to concave 
    Ncnvx = any(qual_<0);
    if any(Ncnvx)
    qual_fix_ = qual_fix;
    qual_fix(:,Ncnvx) = NaN;
    Cc = Cc(Ncnvx(R)); Cr = Cr(Ncnvx(R)); Cl = Cl(Ncnvx(R)); R = R(Ncnvx(R));
    Ifix = qual_(Cc+(R-1)*size(qual_,1))<0; 
    qual_fix(Cl(Ifix)+(R(Ifix)-1)*size(qual_,1)) = qual_fix_(Cl(Ifix)+(R(Ifix)-1)*size(qual_,1));
    qual_fix(Cr(Ifix)+(R(Ifix)-1)*size(qual_,1)) = qual_fix_(Cr(Ifix)+(R(Ifix)-1)*size(qual_,1));
    else
    	allconvex = true;
    end;
    end;
    
    %best candidate element found
    [tmp,CI] = max(qual_fix); CI = CI';
    RI = 1:nbad;
    %give up on removing node, if best candidate element is 
    if (options.consRM || size(edg_,2) == 2) && size(edg_,2) ~= 3 %worse than previously
    Ibetter = qual_(CI+(RI'-1)*size(qual_,1)) > triQtb(badnds); 
    else %inverted (we do not store surface quality in 3D)
    Ibetter = qual_(CI+(RI'-1)*size(qual_,1)) > 0;
    end;
    
    if not(all(Ibetter))
      circles = circles(Ibetter,:); nbad = size(circles,1); 
      qual_   = qual_(:,Ibetter);
      RI = 1:nbad; CI = CI(Ibetter); NN = NN(Ibetter);
      badbadnds = [badbadnds; badnds(not(Ibetter))];
      badnds = badnds(Ibetter);
      badIDs = badIDs(Ibetter);
      if nbad == 0
        break;
      end;
    end;
    Cl = thmp1(CI); I = find(Cl==0); Cl(I) = NN(RI(I));
    Cr = thmp2(CI); I = find(Cr==thmp2(NN(RI))); Cr(I) = 1;
    if size(edg_,2) ~= 0
        edg_ = edg_(Ibetter,:);
    end;
    if size(xy,2) == 2 || size(edg_,2) == 3
    newtri_ = [circles(RI'+(Cr-1)*nbad) circles(RI'+(CI-1)*nbad) circles(RI'+(Cl-1)*nbad)];
    newtriN = [newtriN; badnds];
    newIDs = [newIDs; badIDs];
    qualityN = [qualityN; qual_(CI+(RI'-1)*size(qual_,1))];
    newedg = [newedg; newtri_(NN~=3,[1 3])];
    else %swapping
    qual3Da = qual3Da(:,Ibetter);
    qual3Db = qual3Db(:,Ibetter);
    %newtri_ = [[edg_(:,1) circles(RI'+(Cl-1)*nbad) circles(RI'+(CI-1)*nbad) circles(RI'+(Cr-1)*nbad)]; [edg_(:,2) circles(RI'+(CI-1)*nbad) circles(RI'+(Cl-1)*nbad) circles(RI'+(Cr-1)*nbad)]];
    newtri_ = [[circles(RI'+(Cr-1)*nbad) circles(RI'+(CI-1)*nbad)  circles(RI'+(Cl-1)*nbad) edg_(:,1)]; [circles(RI'+(Cl-1)*nbad) circles(RI'+(CI-1)*nbad)  circles(RI'+(Cr-1)*nbad) edg_(:,2)]];
    newtriN = [newtriN; repmat(badnds,2,1)]; 
    newIDs  = [newIDs;  repmat(badIDs,2,1)];
    qualityN = [qualityN; qual3Da(CI+(RI'-1)*size(qual_,1)); qual3Db(CI+(RI'-1)*size(qual_,1))];
    uppe = [false(nbad,1); true(nbad,1)]; lowe = not(uppe);
    NN1n3 = [NN~=3; false(size(NN))];
    uppe = [false(nbad,1); true(nbad,1)];
    newedg = [newedg; newtri_(NN1n3,[1 3])];
    newfac = [newfac; newtri_(repmat(NN~=3,2,1),[1 3 4]); newtri_(uppe,1:3)];
    newfacN = [newfacN; repmat(badnds(NN~=3),2,1); badnds];
    end;
    newedgN = [newedgN; badnds(NN~=3)];
    newtri = [newtri; newtri_];
    % UPDATE CIRCLES
    circles(RI'+(CI-1)*nbad) = 0;
    qual_(CI+(RI'-1)*size(qual_,1)) = 0;
    qual_(Cr+(RI'-1)*size(qual_,1)) = NaN;
    qual_(Cl+(RI'-1)*size(qual_,1)) = NaN;
    NN = NN - 1;
    I = NN>2; NN = NN(I);
    if not(any(NN>2))
    	break;
    end;
    circles = circles(I,:); qual_ = qual_(:,I); badnds = badnds(I); badIDs = badIDs(I);
    
    if size(edg_,2) == 2
    	qual3Da = qual3Da(:,I);
    	qual3Db = qual3Db(:,I);
    end;
    if size(edg_,2) ~= 0
    	edg_ = edg_(I,:);
    end;
    if size(edg_,2) == 2
    	[circles,qual_,qual3Da,qual3Db] = fix_circles(circles,qual_,qual3Da,qual3Db);
    elseif options.qualRM 
	[circles,qual_] = fix_circles(circles,qual_);
    else
    	qual2 = qual2(:,Ibetter); %if not(all(Ibetter))
    	qual2(CI+(RI'-1)*size(qual2,1)) = 0;
      	qual2(Cr+(RI'-1)*size(qual2,1)) = NaN;
      	qual2(Cl+(RI'-1)*size(qual2,1)) = NaN;
    	qual2 = qual2(:,I);
    	[circles,qual_,qual2] = fix_circles(circles,qual_,qual2);
    end;
end; %while
%end; %ifany
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));
newedg   = newedg( Ikeep(newedgN),:);
newedgN  = newedgN( Ikeep(newedgN));
if size(edg_,2) == 2
newfac   = newfac( Ikeep(newfacN),:);
newfacN  = newfacN(Ikeep(newfacN));
end;

	function [circles,oqual_,oqual2,oqual3] = fix_circles(circles,qual_,qual2,qual3)
ocircles = circles';
[Cg,Rg] = find(ocircles~=0);
[C ,R ] = find(ocircles==0);
move = zeros(size(circles)); move(R+(C-1)*size(circles,1)) = -1; move = cumsum(move')';
Cg = Cg + reshape(move(Rg+(Cg-1)*size(circles,1)),size(Cg));
circles = zeros(size(circles));
circles(Rg+(Cg-1)*size(circles,1)) = ocircles(ocircles~=0);
if nargout >= 2
oqual_ = zeros(size(circles));
oqual_(Rg+(Cg-1)*size(circles,1)) = qual_(find(ocircles)); oqual_ = oqual_';
end;
if nargout >= 3
oqual2 = zeros(size(circles));
oqual2(Rg+(Cg-1)*size(circles,1)) = qual2(find(ocircles)); oqual2 = oqual2';
end;
if nargout >= 4
oqual3 = zeros(size(circles));
oqual3(Rg+(Cg-1)*size(circles,1)) = qual3(find(ocircles)); oqual3 = oqual3';
end;



function [newtri,circles,badnds,badbadnds,newtriN,qualityN,R2R,NN,edg_] = rm_trinds(circles,badnds,badbadnds,xy,Nmetric,triQtb,options,edg_)
NN = sum(circles>0,2);
N3 = NN==3;
R2R = 1:numel(badnds);	
if not(any(N3))
  newtri = []; qualityN = []; newtriN = [];
  return;
end;
if size(xy,2) == 2 || size(edg_,2) == 3
  newtri = circles(N3,[3 2 1]); 
  newtriN = badnds(N3)';
  qualityN = elem_qual(newtri,xy,Nmetric,options); %2D element uninvertible
else %swapping
  newtri = [[edg_(N3,1) circles(N3,1:3)]; ...
            [edg_(N3,2) circles(N3,[2 1 3])]]; 
  qualityN = elem_qual(newtri,xy,Nmetric,options,1);
  newtriN = repmat(badnds(N3)',2,1);
end;
if (options.consRM || size(edg_,2) == 2) && size(edg_,2) ~= 3
  if size(edg_,2) ~= 2
   Ibetter = triQtb(badnds(N3)) < qualityN;
   Ibetter2 = Ibetter;
  elseif size(edg_,2) ~= 3 %swapping
   Ibetter = triQtb(badnds(N3)) < min(reshape(qualityN,[sum(N3) 2]),[],2);
   Ibetter2 = repmat(Ibetter,2,1);
  else %we do not store surface quality in 3D
   Ibetter = 0 < qualityN;
   Ibetter2 = Ibetter;
  end;
  newtri = newtri(Ibetter2,:);
  newtriN = newtriN(Ibetter2);
  qualityN = qualityN(Ibetter2);
  Igood_ = not(N3); Igood_(N3) = Ibetter;
  badbadnds = [badbadnds; badnds(not(Igood_))'];
  badnds = badnds(Igood_); circles = circles(Igood_,:); N3 = N3(Igood_); NN = NN(Igood_);
  if size(xy,2) == 3
  	edg_ = edg_(Igood_,:); edg_ = edg_(not(N3),:); 
  end;
  R2R = (1:numel(badnds))';
end;
circles = circles(not(N3),:); R2R = R2R(not(N3)); NN = NN(not(N3));


function [circles,badnds,edg,badbadnds] = del_dblNcnvx(circles,xy,badnds,edg,options)
NN = sum(circles~=0,2);
[Cc,R] = find(circles'); nbad = size(circles,1);
thmp1 = [0:size(circles,2)-1]'; thmp2 = [2:size(circles,2)+1]';
Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
Cr = thmp2(Cc); I = find(Cr==thmp2(NN(R))); Cr(I) = 1;
I1 = circles(R+(Cr-1)*nbad); I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
I123 = [I1(:) I2(:) I3(:)];
if size(xy,2) == 2
[Ibad,area] = elem_inv(I123,xy);
elseif size(edg,2) == 2 %swapping
[Ibad,area] = elem_inv([[edg(R,1) I123(:,[3 2 1])]; [edg(R,2) I123(:,[2 3 1])] ],xy); 
Ibad = any(reshape(Ibad,size(circles,1),2),2);
else
[Ibad,area] = elem_inv(I123,xy,edg(R,:)); 
end;%Ibad = area < options.minA;
bcircles = false(size(circles));
bcircles(R+(Cc-1)*nbad) = Ibad;
if options.rmnonconvec == 0
	ndblNcnvx = sum(bcircles,2) == 0;
else
	ndblNcnvx = sum(bcircles,2) < 2;
end;
circles = circles(ndblNcnvx,:);
badbadnds = badnds(not(ndblNcnvx));
badnds = badnds(ndblNcnvx);
if size(xy,2) == 3
	edg    = edg(ndblNcnvx,:);
end;

function [newtri,qualityN,badnds,newtriN,badbadnds,newIDs] = prag_coarse(circles,Nmetric,xy,badnds,badIDs,triQtb,NN,options,nvec)
badbadnds = zeros(0,1); badnds_start = badnds;
newtri = zeros(0,3); newtriN = []; qualityN = []; newIDs = [];

[Cc,R] = find(circles'~=0); nbad = size(circles,1);
edg = [badnds(R) reshape(circles(R+(Cc-1)*nbad),size(R))];
edgL = inf(size(circles));
edgL(R+(Cc-1)*nbad) = elem_qual(edg,xy,Nmetric,options);
[edgL_,Iedg] = sort(edgL,2);
Iedg(isinf(edgL)) = 0;
theC = 1;

while true
[Cc,R] = find(circles'~=0); nbad = size(circles,1);
thmp1 = [0:size(circles,2)-1]'; Cl = thmp1(Cc); I = find(Cl==0); Cl(I) = NN(R(I));
I2 = circles(R+(Cc-1)*nbad); I3 = circles(R+(Cl-1)*nbad);
C_ = Iedg(R+(theC-1)*nbad);
Ishrt = circles(R+(C_(:)-1)*nbad);
newtri_ = [Ishrt(:) I2(:) I3(:)];
[newtri__,tmp] = sort(newtri_,2);
Igood = and(newtri__(:,1) ~= newtri__(:,2),newtri__(:,2) ~= newtri__(:,3));
newtri_ = newtri_(Igood,:); R = R(Igood); Cc = Cc(Igood);
if size(xy,2) == 2
   newqual = elem_qual(newtri_,xy,Nmetric,options,1);
else
   newqual = elem_qual(newtri_,xy,Nmetric,options,nvec(R,:));
end;
if options.consRM && size(nvec,2) ~= 3
	Igood = newqual > triQtb(badnds(R));
else
	Igood = newqual > options.minqual;
end;
allgood = true(size(circles)); allgood(R+(Cc-1)*nbad) = Igood;
Igood = all(allgood,2); Ibad = not(Igood);
newtri = [newtri; newtri_(Igood(R),:)];
qualityN = [qualityN; newqual(Igood(R))];
newtriN  = [newtriN ; badnds(R(Igood(R)))];
newIDs   = [newIDs  ; badIDs(R(Igood(R)))];
%update circles
circles = circles(Ibad,:); badnds = badnds(Ibad); NN = NN(Ibad); Iedg = Iedg(Ibad,:); badIDs = badIDs(Ibad);
if size(xy,2) == 3
	nvec = nvec(Ibad,:);
end; 
if size(circles,1) == 0
	break;
end;
theC = theC+1;
Ibad = NN < theC; Igood = not(Ibad);
badbadnds = [badbadnds; badnds(Ibad)];
circles = circles(Igood,:); badnds = badnds(Igood); NN = NN(Igood); Iedg = Iedg(Igood,:); badIDs = badIDs(Igood);
if size(xy,2) == 3
	nvec = nvec(Igood,:);
end;
if size(circles,1) == 0
	break;
end; %size(circles)
end; %badbadnds
Ikeep = true(max(badnds_start),1); Ikeep(badbadnds) = false;
badnds = badnds_start(Ikeep(badnds_start));
%next four lines are irrelevant outside the context of region IDs
newtri   = newtri(  Ikeep(newtriN),:);
qualityN = qualityN(Ikeep(newtriN));
newIDs = newIDs(Ikeep(newtriN));
newtriN  = newtriN( Ikeep(newtriN));


%hold on; trimesh(newtri,xy(:,1),xy(:,2)); hold off;