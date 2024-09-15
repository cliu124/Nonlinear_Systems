function fb=bgradu2f(p,f,b,u)
% BGRADU2F: include "b grad u" into f for error-est and mesh-adapt.
%
%  fb=bgradu2f(p,f,b,u)
% f,b as returned from p.f 
%
% See also pdegrad, meshada, errcheck, meshref
[ux,uy]=pdegrad(p.mesh.p,p.mesh.t,p.mat.fill*u(1:p.nu));
fb=f;n=p.nc.neq; 
for i=1:n
    for j=1:n
        fb(i,:)=fb(i,:)+b(2*n*(j-1)+2*(i-1)+1,:).*ux(j,:)...
            +b(2*n*(j-1)+2*(i-1)+2,:).*uy(j,:);
    end
end
