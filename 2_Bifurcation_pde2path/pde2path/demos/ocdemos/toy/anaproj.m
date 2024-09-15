function [eigvxy,muv,F,F2]=anaproj(p,u)
% compute the monodromy matrix and projections for the ODE toy system 
% via polar coordinates in the states
par=p.u(p.nu+1:end); syms x y r y1;
% derivatives of polar coordinates
dphi=@(x0,y0) subs([diff(atan(y/x),x);diff(atan(y/x),y)],{x,y},{x0,y0});
dr=@(x0,y0) subs([diff(sqrt(x^2+y^2),x);diff(sqrt(x^2+y^2),y)],{x,y},{x0,y0});
% compute eigenvalues in polar coordinates at each time step
[a,b]=eig([par(3)*(-1+3*y1*r^2) 0 par(3)*r^3 0;0 0 0 0;0 0 0 1*par(1);0 0 2*pi*cos(2*pi*y1)*par(1) 0]);
a=simplify(a,10);
b=subs(exp(2*pi*b),{r,y1},{sqrt(u(1)^2+u(2)^2),u(3)}); muv=double(diag(b)); % floquet-multipliers
eigvxy=zeros(4,4); % compute eigenvectors in cartesian coordinates
for k=1:4
    eigv=subs(a(:,k),{r,y1},{sqrt(u(1)^2+u(2)^2),u(3)});
    eigv=simplify(eigv,10);
    if double(eigv(1))~=0
        eigvxy(1:2,k)=eigv(1)*dr(u(1),u(2));
    end
    if double(eigv(2))~=0
        eigvxy(1:2,k)=eigv(2)*dphi(u(1),u(2));
    end
    eigvxy(3:4,k)=eigv(3:4);
end
eigvxy=double(eigvxy);
unst=eigvxy(:,muv<1-1e-6);
F2=null(unst.').'; % projection on stable eigenspace
unst=eigvxy(:,muv<1+1e-6); % projection on center-unstable eigenspace
F=null(unst.').';
