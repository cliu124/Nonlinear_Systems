function out=cmcbra(p,u)
par=u(p.nu+1:end); X=p.X;  %del=par(10)/p.pold;
%a=[ones(p.np,2) del*ones(p.np,1)]; X=a.*p.X;
V=getV(p,u); A=getA(p,u);
q=boundary_faces(p.tri); ids=unique([q(:,1);q(:,2)]); 
x1=min(X(ids,1)); x2=max(X(ids,1)); 
y1=min(X(ids,2)); y2=max(X(ids,2)); 
z1=min(X(ids,3)); z2=max(X(ids,3)); 
r=min(sqrt(X(:,1).^2+X(:,2).^2)); 
R=max(sqrt(X(:,1).^2+X(:,2).^2)); 
out=[par; V;A;x1;x2;y1;y2;z1;z2;r;R]; 