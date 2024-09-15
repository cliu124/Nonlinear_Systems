function [crnds,edg,edg2,edg2tri,edg2ID] = geom_crnds(bndmesh,addcrs)
if not(isfield(bndmesh,'fac'))
nds = bndmesh.edg(:);
IDs = repmat(bndmesh.IDs,2,1);
[nds,I] = sort(nds);
IDs = IDs(I);
I2 = [false; and(nds(1:end-1)==nds(2:end),IDs(1:end-1)~=IDs(2:end))];
crnds = nds(I2);
%I3_ = find(I3);
%snglnds = nds(I3_(diff([0; I3_]) == 1));
snglnds = nds([nds(1) ~= nds(2); and(nds(1:end-2)~=nds(2:end-1),nds(2:end-1)~=nds(3:end)); nds(end-1) ~= nds(end)]);
crnds = [crnds; snglnds]; %an edge can stop in 2D space, but we do not check for isolated nodes
%bndmesh.crnds = [bndmesh.crnds; crnds];
else %3D
nds = bndmesh.fac(:);
IDs = repmat(bndmesh.IDs,3,1);
%get corner nodes
[edg2,edg2ID,edg2tri] = books(bndmesh.fac,bndmesh.IDs);
I = edg2ID(:,1)~=edg2ID(:,2);
if size(edg2ID,2) > 2
I = or(I,sum(edg2ID~=0,2)>2);
end;
edg = edg2(I,:);
if nargin == 2
%edg2ID = edg2ID(I,:);
%edg2tri = edg2tri(I,:);
edgs = [1:size(edg,1) 1:size(edg,1)]';
[nds,I] = sort(edg(:));
%IDs = IDs(I);
edgs = edgs(I);
I2 = [false; and(and(nds(2:end-1)==nds(3:end),edgs(2:end-1)~=edgs(3:end)),and(nds(2:end-1)==nds(1:end-2),edgs(2:end-1)~=edgs(1:end-2))); false];
crnds = nds(I2);
if numel(crnds) ~=0
crnds = crnds([crnds(1:end-1)~=crnds(2:end); true]);
%bndmesh.crnds = [bndmesh.crnds; crnds];
else
crnds = [];
end;
else
crnds = [];
end;
%bndmesh.edg = edg;
%bndnds2 = false(max(edg(:)),1); bndnds2(edg(:)) = true; bndnds2 = find(bndnds2);
end;

 %a face can stop in 3D space, but we do not check for isolated nodes or lines going through 3D space

function [edg,edg2ID,edg2tri] = books(tria,IDs)
edga = [tria(:,2) tria(:,1); tria(:,3) tria(:,2); tria(:,1) tria(:,3)];
edga2tri = [1:size(tria,1) 1:size(tria,1) 1:size(tria,1)]';
IDs = repmat(IDs,3,1);
edg = sort(edga,2); %biggest in last column
[edg,Is] = sortrows(edg);
tris = edga2tri(Is);
IDs = IDs(Is);
d =[true; or(edg(1:end-1,1)~=edg(2:end,1), ...
            edg(1:end-1,2)~=edg(2:end,2))]; Nd = sum(d);
edg = edg(d,:);
edgs = cumsum(d);
edg2tri = rpval2M(edgs,tris); edg2tri = rpval2M_clean(edg2tri);
edg2ID = rpval2M(edgs,IDs); edg2tri = rpval2M_clean(edg2tri);