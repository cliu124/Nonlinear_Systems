function bndmesh = bndmesh_polygon(tri,xy,bndmesh,options)
[edg,edg2tri,tri2edg,nd2tri,nd2edg] = bks_all(tri);
edg = edg(edg2tri(:,2)==0,:);
x = xy(:,1); y = xy(:,2);
c = [mean(x(edg),2) mean(y(edg),2)];
t = xy(edg(:,1),:)-xy(edg(:,2),:);
t = t./repmat(sqrt(t(:,1).^2+t(:,2).^2),1,2);
IDs = zeros(size(edg,1),1);
while true
    [C,n] = min(IDs);
    IDs(n) = max(IDs) + 1;
    I = 1:numel(IDs);
    notnset = I ~= n;
    dists = abs(t(notnset,2).*(xy(edg(notnset,1),1) - xy(edg(n,1),1))-...
                t(notnset,1).*(xy(edg(notnset,1),2) - xy(edg(n,1),2))) < options.geomtol;
    angles = 1-abs(t(notnset,1).*t(n,1)+t(notnset,2).*t(n,2))<options.geomtol; 
    % angles = arccos(abs(t[notnset,0]*t[n,0]+t[notnset,1]*t[n,1]))<1e-12
    I_ = I(notnset);
    IDs(I_(and(angles,dists))) = IDs(n);
    if all(IDs ~= 0)
        break;
    end;
end;
if options.verbose
 disp(sprintf('Found %i co-linear edges', max(IDs)));
end;
bndmesh.edg = edg;
bndmesh.IDs = IDs;
if isfield(bndmesh,'crnds')
bndmesh.crnds = [bndmesh.crnds; geom_crnds(bndmesh)];
else
bndmesh.crnds = geom_crnds(bndmesh);
end;
