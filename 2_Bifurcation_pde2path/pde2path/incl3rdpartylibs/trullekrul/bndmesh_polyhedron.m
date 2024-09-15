function bndmesh=bndmesh_polyhedron(tri,xy,bndmesh,options)
[fac,fac2tri,tri2fac]=bks_all3D(tri);
fac=fac(fac2tri(:,2)==0,:);
x=xy(:,1); y=xy(:,2); z=xy(:,3);
c=[mean(x(fac),2) mean(y(fac),2) mean(z(fac),2)];
v1=xy(fac(:,1),:)-xy(fac(:,2),:);
v2=xy(fac(:,3),:)-xy(fac(:,2),:);
n=cross(v1,v2,2); %[v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) ...
     %v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) ...
     %v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)];
n=n./repmat(sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2),1,3);
IDs=zeros(size(fac,1),1);
while true    
    [C,nn]=min(IDs); 
    IDs(nn)=max(IDs) + 1;
    I=1:numel(IDs);
    notnset=I ~= nn;
    dists=abs(n(notnset,1) .* (xy(fac(notnset,1),1) - xy(fac(nn,1),1)) + ...
      n(notnset,2) .* (xy(fac(notnset,1),2) - xy(fac(nn,1),2)) + ...
      n(notnset,3) .* (xy(fac(notnset,1),3) - xy(fac(nn,1),3))) < options.geomtol;
    angles=1-abs(n(notnset,1).*n(nn,1) + ...
      n(notnset,2).*n(nn,2) + n(notnset,3).*n(nn,3))<options.geomtol; 
    I_=I(notnset);
    IDs(I_(and(angles,dists)))=IDs(nn);
    if all(IDs ~= 0)
        break;
    end;
end;
if options.verbose
 disp(sprintf('Found %i co-linear faces', max(IDs)));
end;
bndmesh.fac=fac;
bndmesh.IDs=IDs;
if isfield(bndmesh,'crnds')
bndmesh.crnds=[bndmesh.crnds; geom_crnds(bndmesh,1)];
else
bndmesh.crnds=geom_crnds(bndmesh,1);
end;
