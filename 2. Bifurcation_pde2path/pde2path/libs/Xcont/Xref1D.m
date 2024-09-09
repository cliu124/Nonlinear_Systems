function [p,un,taun]=Xref1D(p,idx)
ne=length(idx); ids=sort(idx); uo=p.u(1:p.nu); tauo=p.tau(1:p.nu); Xo=p.X; onp=p.np; 
for i=1:ne % loop over edges to refine 
    j=ids(i); i1=p.tri(j,1); i2=p.tri(j,2); 
    i1=i1+i-1; if i2~=1; i2=i2+i-1; end % correct for already inserted points
    x1=p.X(i1,:); x2=p.X(i2,:); xm=0.5*(x1+x2); % new point 
    if i2~=1; p.X=[p.X(1:i1,:); xm; p.X(i2:end,:)]; % insert, if not last elem 
    else p.X=[p.X(1:i1,:); xm]; % otherwise just append
    end 
end
% interpolate u, tau 
xo=Xo(:,1); yo=Xo(:,2); 
xn=p.X(:,1); yn=p.X(:,2); np=size(p.X,1); p.np=np; p.nt=np; 
p.tri(1:np,1)=1:np; p.tri(1:np,2)=circshift(1:np,-1); 
for i=1:p.nc.neq % interpolate u and tau to new mesh 
    un((i-1)*np+1:i*np)=p2interpol(xn,yn,uo((i-1)*onp+1:i*onp),xo,yo,p); 
    taun((i-1)*np+1:i*np)=p2interpol(xn,yn,tauo((i-1)*onp+1:i*onp),xo,yo,p); 
    upn((i-1)*np+1:i*np)=p2interpol(xn,yn,p.up((i-1)*onp+1:i*onp),xo,yo,p); 
end
p.up=upn'; un=un'; taun=taun'; 
return 
size(uo), size(un)
xo(1:12)'
xn(1:12)'
uo(onp+1:onp+12)'
un(np+1:np+12)'