function [edg2clr,clr2edg] = bks_clr(edg_ngh,options) 
%%% ALGORITHMS
%options.greedyCLR == 3  greedy with Largest degree first post-processing
%options.greedyCLR == 2  greedy with Jones Plassmann post-processing
%options.greedyCLR == 1  greedy (colouring)
%options.greedyCLR == 0  Jones Plassmann
%options.greedyCLR == -1 Largest degree first 

if options.debug == 2 %check input
for i=1:size(edg_ngh,1)
if not(all(any(edg_ngh(edg_ngh(i,edg_ngh(i,:)~=0),:) == i,2)))
	error(sprintf('input to coloring algorithm is flawed (%0.0f)',i));
end; end; end;

edgs = size(edg_ngh,1);
if not(options.greedyCLR>0) && mean(sum(edg_ngh~=0,2)==size(edg_ngh,2)) < 0.5
	edg2clr = clr_JP(edg_ngh,options);
else
mxtrs = 70;
nghI = edg_ngh~=0;
trs = 0; clrs = size(edg_ngh,2)+(mean(sum(edg_ngh~=0,2)==size(edg_ngh,2)) > 0.5);
edg2clr = ceil(rand(edgs,1)*clrs); %1+mod(1:edgs,clrs);
while 1
    if trs > mxtrs
        trs = 0;
        clrs = clrs + 1;
        warning(sprintf('### EDGE GROUPING: GROUP SIZE INCREASED(%0.0d) ###',clrs-1));
    end;
    if trs == 0
        badedg = (1:edgs)';
    end;
    badedgclr = edg_ngh(badedg,:); badedgclr(nghI(badedg,:)) = edg2clr(badedgclr(nghI(badedg,:)));
    cmpclr = edg2clr(repmat(badedg,1,size(edg_ngh,2)));
    I = sum(cmpclr == badedgclr,2)>0; badedg = badedg(I);
    if numel(badedg) == 0
        break;
    end;
    % SUGGEST NEW COLORS TO BAD EDGES
    trs = trs + 1;
    I = find(cmpclr(I,:) == badedgclr(I,:));
    badedg_ = 1:edgs+1; badedg_(badedg) = 0; 
    tknclrs = edg_ngh(badedg,:); tknclrs(nghI(badedg,:)) = badedg_(tknclrs(nghI(badedg,:)));
    
    %do not consider colours of badedgs as taken
    tknclrs(tknclrs~=0) = edg2clr(tknclrs(tknclrs~=0));
    if max(sum(tknclrs ~= 0,2)) == clrs
        trs = mxtrs + 1;
        continue;
    end;
    Nbad = numel(badedg);
    llclrs = repmat((1:clrs),Nbad,1); 
    [R,tmp,C] = find(tknclrs); 
    llclrs(R+(C-1)*Nbad) = 0;
    vlclrs = sum(llclrs~=0,2); llclrs = sort(llclrs,2);
    C = clrs+1-ceil(rand(Nbad,1).*vlclrs)'; nwclrs = llclrs((1:Nbad)+(C-1)*Nbad);
    edg2clr(badedg) = nwclrs;
end;

if options.greedyCLR > 1
%REDUCE CROMATIC NUMBER
if options.greedyCLR == 3
Nn = sum(edg_ngh~=0,2);
end;

while true
done = 1;
for i=1:clrs
	inds = edg2clr==i; %inds = edg2clr_(1:end-1)==i;
	finds = find(inds); 
	Ni = numel(finds);
	if Ni == 0
		continue;
	end;
	tknclrs = edg_ngh(inds,:); tknclrs(nghI(inds,:)) = edg2clr(tknclrs(nghI(inds,:))); 
	tknclrs = [edg2clr(finds) tknclrs];
	llclrs =  repmat((1:clrs),Ni,1); 
	[R,C,Cc] =  find(tknclrs);
	llclrs(R+(Cc-1)*Ni) = clrs+1;
	minclr = min(llclrs,[],2);
	if options.greedyCLR == 3
	%COLOUR MOST NEIGHBOUGHS FIRST
	llclrs =  repmat((1:clrs),Ni,1); 
	llclrs(R+(Cc-1)*Ni) = 0;
	maxclr = max(llclrs,[],2);
	nNngh = edg_ngh(inds,:); nNngh(nghI(inds,:)) = Nn(nNngh(nghI(inds,:)));
	nNngh = reshape(nNngh,Ni,size(edg_ngh,2));  %reshaping when Ni==1
	I2 = repmat(Nn(finds),1,size(nNngh,2)) < nNngh;
	bgngh = zeros(size(nNngh)); bgngh(I2) = tknclrs(I2);
	I = and(and(max(bgngh,[],2)<minclr,minclr < edg2clr(finds)),Nn(finds)~=0);
	Im = and(Nn(finds)~=0,not(I));
	edg2clr(finds(Im)) = maxclr(Im);
	Nn(finds(I)) = 0;
	else
	I = minclr < edg2clr(finds);
	end;
	edg2clr(finds(I)) = minclr(I);	
	if any(I)
		done = 0;
	end;
end; %for
if done == 1
	break;
end;
end;%while true
end; %reduce cromatic number
end; %greedyCLR


if nargout == 2
    clr2edg = env_table(edg2clr);
end;
if options.debug
    % SIMPLE CHECK
    edg2clr_ = NaN(size(edg_ngh)); edg2clr_(edg_ngh~=0) = edg2clr(edg_ngh(edg_ngh~=0));
    cmpmat = repmat(edg2clr(1:edgs),1,size(edg_ngh,2)); cmpmat = edg2clr_ == cmpmat;
    if any(cmpmat(:)) %0~=sum(sum(cmpmat))
        [R,C] = find(cmpmat); 
        error(sprintf('EGDE NUMBER %0.0f HAS WRONG COLOUR',R(1)))
    end;
end;


function edg2clr = clr_JP(edg_ngh,options)
edgs = size(edg_ngh,1);
nghI = edg_ngh~=0;
nN = sum(edg_ngh~=0,2);
clrs = size(edg_ngh,2)+1;
edg2clr = zeros(edgs,1);
badedg = (1:edgs)';
while numel(badedg) ~= 0
    badedgf_ = [zeros(edgs,1)];
    badedgf_(badedg) = badedg;
    %edges to color (they have the largest number):
    badedg_ngh = edg_ngh(badedg,:); badedg_ngh(nghI(badedg,:)) = badedgf_(badedg_ngh(nghI(badedg,:)));
    if options.greedyCLR == 0
	    colorI = badedg > max(badedg_ngh,[],2); 
    else
	    %edges to color have the most neighbours (largest number in case of conflicts):
	    badedgN = badedgf_; badedgN(badedg) = nN(badedg);
	    badedg_nghN = edg_ngh(badedg,:); badedg_nghN(nghI(badedg,:)) = badedgN(badedg_nghN(nghI(badedg,:)));
	    [maxN,C] = max(badedg_nghN,[],2); I = repmat(maxN,1,size(edg_ngh,2)) == badedg_nghN; 
	    maxNedg = zeros(size(badedg_ngh)); maxNedg(I) = badedg_ngh(I); maxNedg = max(maxNedg,[],2);
	    colorI = or(nN(badedg) > maxN,and(nN(badedg) == maxN,badedg>maxNedg));
    end;
        
    badedgI = badedg(colorI); 
    badedg_ = 1:edgs; badedg_(badedgI) = 0; 
    tknclrs = edg_ngh(badedgI,:); tknclrs(nghI(badedgI,:)) = badedg_(tknclrs(nghI(badedgI,:)));
    %do not consider colours of badedgs as taken
    tknclrs(tknclrs~=0) = edg2clr(tknclrs(tknclrs~=0));

    tknclrs = [zeros(numel(badedgI),clrs-size(tknclrs,2)) tknclrs];
    R = repmat((1:numel(badedgI))',1,size(tknclrs,2)); I = find(tknclrs); C = tknclrs(I);
    tknclrs = zeros(size(tknclrs)); tknclrs(R(I)+(C-1)*numel(badedgI)) = C;
    llclrs = repmat((1:clrs),numel(badedgI),1); llclrs(find(tknclrs)) = 0;
    vlclrs = sum(llclrs>0,2); llclrs = sort(llclrs,2);
    C = clrs-vlclrs'+1; nwclrs = llclrs((1:numel(badedgI))+(C-1)*numel(badedgI));
    edg2clr(badedgI) = nwclrs;
    clrs = max(clrs,max(nwclrs)+1);
    badedgf_(badedgI) = 0; badedg = badedgf_(find(badedgf_));
end;