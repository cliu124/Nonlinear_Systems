function [p,un,taun]=Xcoarse1D(p,idx)
% Xcoarse1D: coarsen X (curve) by deleting points  
uo=p.u(1:p.nu); tauo=p.tau(1:p.nu); Xo=p.X; onp=p.np; 
nl=setdiff(1:onp,idx); p.X=p.X(nl,:); % delete points 
xo=Xo(:,1); yo=Xo(:,2); 
xn=p.X(:,1); yn=p.X(:,2); np=size(p.X,1); p.np=np; p.nt=np; 
p.tri=[]; p.tri(1:np,1)=1:np; p.tri(1:np,2)=circshift(1:np,-1); 
for i=1:p.nc.neq % interpolate u and tau to new mesh 
    un((i-1)*np+1:i*np)=p2interpol(xn,yn,uo((i-1)*onp+1:i*onp),xo,yo,p); 
    taun((i-1)*np+1:i*np)=p2interpol(xn,yn,tauo((i-1)*onp+1:i*onp),xo,yo,p); 
    upn((i-1)*np+1:i*np)=p2interpol(xn,yn,p.up((i-1)*onp+1:i*onp),xo,yo,p); 
end
p.up=upn'; un=un'; taun=taun'; 