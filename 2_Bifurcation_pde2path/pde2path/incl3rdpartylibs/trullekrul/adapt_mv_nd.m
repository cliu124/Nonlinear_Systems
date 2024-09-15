function [xy,triQ,bks,ndone] = adapt_mv_nd(xy,tri,Nmetric,bndmesh,triQ,bks,geomfunc,options,repeat)
ndone = 0;
nds = size(xy,1);
if nargin == 8
	repeat = 1;
end;
if size(xy,2) == 3
	[tmp,edg1,edg2] = geom_crnds(bndmesh); bndmesh.edg = edg1;
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	bndnds2_ = false(size(xy,1),1); bndnds2_(bndmesh.fac(:)) = true;
	if size(Nmetric,2) ~= 6
	Nmetric = Nmetric(:,1:6);
	end;
else
	bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.edg(:)) = true;
	if size(Nmetric,2) ~= 3
	Nmetric = Nmetric(:,1:3);
	end;
end;

nd2tri = inv_table(tri);
ngh = bks_nd2tri2nd2(tri,nd2tri);

nd2clr = bks_clr(ngh,options); nd2clr(bndmesh.crnds) = 0;
if isfield(options,'TS')
	if size(xy,2) == 3
		nd2clr(bndmesh.fac(bndmesh.IDs==options.TS,:)) = 0;
	else
		nd2clr(bndmesh.edg(bndmesh.IDs==options.TS,:)) = 0;
	end;
end;
clrs = max(nd2clr);


%BOUNDARY NODES ARE AVERAGED WITH NEIGHBORING BOUNDARY NODES ONLY
if size(xy,2) == 3
bndngh = gen_bndngh3(bndmesh,bndnds2_); 
ngh(bndnds2_,1:size(bndngh,2)) = bndngh;
ngh(bndnds2_,size(bndngh,2)+1:end) = 0;
end;

if any(bndnds_)
bndngh = gen_bndngh2(bndmesh,bndnds_); ngh(bndnds_,[1 2]) = bndngh;
ngh(bndnds_,3:end) = 0;
end;
%bndngh = gen_bndngh2(bndmesh,nds);  ngh(bndnds_(1:end-1),[1 2]) = bndngh(bndnds_(1:end-1),:);

if size(xy,2) == 3
	bndnds1_ = bndnds_; bndnds_ = bndnds2_;
	bndmesh = rmfield(bndmesh,'edg');
end;

NN = sum(ngh~=0,2); 
ngh_ = ngh;
%start moving
for jj=1:repeat
for i=1:clrs
    inds = nd2clr==i;
    if options.mntn_bks
    	inds = and(inds,bks.mvnd);
    end;
    posnewi = calc_posnew(xy,ngh(inds,:),NN(inds),Nmetric,inds,options);
    if options.MVspd ~= 1.
    	posnewi = posnewi*options.MVspd+xy(inds,:)*(1-options.MVspd);
    end;
    posnew = xy; posnew(inds,:) = posnewi;
    % FIX CURVED BOUNDARY 
    indsB = and(inds,bndnds_);
    if any(indsB) && numel(geomfunc) ~= 0
    	  if size(xy,2) == 2
            dist = geomfunc(posnew(indsB,:));
            else
            indsB1 = and(indsB,not(bndnds1_));
            indsB2 = and(indsB,bndnds1_);
	  dist = zeros(sum(indsB),4);
	  dist(indsB1(indsB),:) = geomfunc{1}(posnew(indsB1,:));
	  if any(bndnds1_)
	  dist(indsB2(indsB),:) = geomfunc{2}(posnew(indsB2,:));
	  end;
            end;
            posnew(indsB,:) = posnew(indsB,:) - dist(:,2:size(xy,2)+1) .* repmat(dist(:,1),1,size(xy,2));	
    end;
    if numel(geomfunc) ~= 0
    	%we will not allow interior nodes outside the geometry:
    	indsnB = and(inds,not(bndnds_));
    	if size(xy,2) ~= 2
    	  dist = geomfunc{1}(posnew(indsnB,:));
    	else
    	  dist = geomfunc(posnew(indsnB,:));
    	end;
    	good = dist(:,1) > -10.*eps;
    	if any(not(good))
    	warning(sprintf('SMOOTHING CONSTRAINED DUE TO ATTEMPTED CREATION OF %0.0f INTERIOR NODE(S) OUTSIDE GEOMETRY',sum(not(good))));
    	inds(indsnB) = good;
    	%moveback = false(size(xy,1),1); moveback(indsnB) = not(good);
    	%posnew(moveback,:) = xy(moveback,:);
    	end;
    end;

    tomove = inds;
    
    qualold  = gen_maxedgE(xy   ,ngh(tomove,:),tri,nd2tri(tomove,:),Nmetric,bndnds_,triQ,options,1); %old element qualities

    qualnew = gen_maxedgE(posnew,ngh(tomove,:),tri,nd2tri(tomove,:),Nmetric,bndnds_,triQ,options,0); %new element qualities
    if options.consMV
	worstimproved = min(qualnew,[],2) > min(qualold,[],2);
    else
    	worstimproved = min(qualnew,[],2) > options.minqual; %accepting all (valid) smoothing
    end;
    %update coordinates and triQ
    tomove(tomove) = worstimproved;
    if sum(tomove) ~= 0
        ndone = ndone + sum(tomove);
    end;
    xy(tomove,:) = posnew(tomove,:);
    if options.mntn_bks
    bks.mvnd(and(inds,not(tomove))) = false;
    bks.rmnd(tomove) = true; %mvnd is already true
    chnd = ngh(tomove,:); chnd = chnd(chnd~=0);
    bks.rmnd(chnd) = true;
    bks.mvnd(chnd) = true;
    end;
    
    qualnew = qualnew(worstimproved,:);
    afftri  = nd2tri(tomove,:);
    I = afftri ~= 0;
    triQ(afftri(I)) = qualnew(I);
    if options.debug
	sanity_check(tri,xy,triQ,Nmetric,options);
	if numel(geomfunc) ~= 0 && size(xy,2) == 3
bndnds_ = false(size(xy,1),1); bndnds_(bndmesh.fac) = true;
dist = geomfunc(xy(bndnds_,:));
if max(abs(dist(:,1))) > options.geomtol
	error('boundary node far from boundary');
end;
end;
    end;
end;
end;
ndone = ndone / size(xy,1);

function [qual,afftri] = gen_maxedgE(xye,ngh,tri,nd2tri,Nmetric,bndnds_,triQ,options,triQuse)
%qual is matrix (nds,1), that gives the qualities for elements with the nodes 
    I = nd2tri~=0; [R,C] = find(I); nds = sum(I(:));
    afftri_ = nd2tri(R+(C-1)*size(nd2tri,1));
    if triQuse
	    quality = triQ(afftri_);    	
    else
	    ntri = [reshape(tri(afftri_,:),nds,size(tri,2))];  %reshape, when nds==1?
	    quality = elem_qual(ntri,xye,Nmetric,options,1);
    end;
qual   = Inf(size(nd2tri));         qual(R+(C-1)*size(qual,1)) = quality;
if nargout == 2
    afftri = zeros(size(nd2tri)); afftri(R+(C-1)*size(qual,1)) = afftri_;
end;


function posnewi = calc_posnew(xy_,ngh_,NN_,Nmetric,inds,options)
dngh = size(ngh_);
I = ngh_~=0; 
nghI = ngh_(I);
[R,C] = find(I);
self = repmat(reshape(find(inds)',dngh(1),1),1,dngh(2));
edg_ = [reshape(self(I),sum(I(:)),1) reshape(ngh_(I),sum(I(:)),1)];
L = elem_qual(edg_,xy_,Nmetric,options);
Lo = zeros(dngh); Lo(I) = L;
x = zeros(size(ngh_)); x(I) = xy_(nghI,1);
y = zeros(size(ngh_)); y(I) = xy_(nghI,2);
posnewi = [sum(Lo.*x,2) sum(Lo.*y,2)]./repmat(sum(Lo,2),1,2);
if size(xy_,2) == 3
z = zeros(size(ngh_)); z(I) = xy_(nghI,3);
posnewi = [posnewi sum(Lo.*z,2)./sum(Lo,2)];
end;



function ngh = gen_bndngh2(bndmesh,bndnds_)
nbnd = sum(bndnds_);
bndnds = zeros(size(bndnds_,1),1); bndnds(bndnds_) = 1:nbnd; 
bndndsI = [find(bndnds_); 0];
nd2edg = inv_table(bndnds(bndmesh.edg));
nd2edg = nd2edg(:,[1 2]); %we will not move corners anyway
edg = [bndnds(bndmesh.edg); [0 NaN]]; nd2edg(nd2edg==0) = size(edg,1); %lone nodes only have a single edge
self = repmat((1:nbnd),2,1);
edg1 = edg(nd2edg(:,1),:)';
edg2 = edg(nd2edg(:,2),:)';
ngh = [edg1(edg1 ~=  self) edg2(and(edg2 ~= self, not(isnan(edg2))))]; 
ngh(ngh==0) = nbnd+1;
ngh = bndndsI(ngh);

function ngh = gen_bndngh3(bndmesh,bndnds_);
nbnd = sum(bndnds_);
bndnds = zeros(size(bndnds_,1),1); bndnds(bndnds_) = 1:nbnd; 
bndndsI = [find(bndnds_); 0];
nd2fac = inv_table(bndnds(bndmesh.fac));
ngh = bks_nd2tri2nd2(bndnds(bndmesh.fac),nd2fac);
ngh(ngh==0) = nbnd+1;
ngh = bndndsI(ngh);