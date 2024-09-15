function [dlam,al]=hogetnf(s,A,om,mu0,q,del) 
% HOGETNF: get dlam and alpha for Hopf initial guess  
%
% [dlam,al]=hogetnf(s,A,om,mu0,q,del) 
% from the NF coefficients for Hopf bifurcation, following [Kuz04, §10.2.2] 
% s=standard p2p-struct, avoiding p since is that used for adjoint Evec
% A=G_u, (om,mu0)=imag/real part at bif, q=Evec, del for FD
M=getM(s); % om=imag(muv(1)); mu0=real(muv(1)); 
u=s.u; dellam=1e-4; del2=del^2; % q=q/norm(q,2); 
% step 1
qr=real(q); qi=imag(q); j0=s.sol.j0; % qr'*qr+qi'*qi, qr'*qi,
s.nc.eigref(j0)=-om*1i; s.nc.neig(j0)=1; % compute p at -om*i
[ineg,muv,Va]=spcalc(A',s,j0);% om, j0, imag(muv(1)), pause
p=Va(:,1); pr=real(p); pi=imag(p); 
npar=size(u,1)-s.nu; qz=zeros(npar,1); % for up fill-up of par 
th=angle(pr'*qr+pi'*qi+(pr'*qi-pi'*qr)*1i); 
p=p*exp(1i*th); pr=real(p); pi=imag(p); 
pnorm=pr'*qr+pi'*qi; pr=pr/pnorm; pi=pi/pnorm; 
%qr'*qr+qi'*qi, qr'*qi, pr'*qr+pi'*qi, pr'*qi-pi'*qr % check normalization 
%norm(A*qr+om*M*qi), norm(-om*M*qr+A*qi), 
%norm(A'*pr-om*M*pi), norm(om*M*pr+A'*pi), pause
% step2
f0=nodalf(s,u); fp=nodalf(s,u+del*[qr;qz]); 
fm=nodalf(s,u-del*[qr;qz]); 
a=(fm-2*f0+fp)/del2; %norm(a)
fp=nodalf(s,u+del*[qi;qz]); fm=nodalf(s,u-del*[qi;qz]); 
b=(fm-2*f0+fp)/del2; %norm(b)
f1p=nodalf(s,u+del*[qr+qi;qz]); f2p=nodalf(s,u+del*[qr-qi;qz]);
f1m=nodalf(s,u-del*[qr+qi;qz]); f2m=nodalf(s,u-del*[qr-qi;qz]);
c=0.25*(f1m+f1p-f2m-f2p); %norm(c), pause
% step 3
r=A\(M*(a+b)); sv=(-A+2i*M*om)\(M*(a-b+2i*c)); 
sr=real(sv); si=imag(sv); 
% step 4
f1p=nodalf(s,u+del*[qr+r;qz]); f2p=nodalf(s,u+del*[qr-r;qz]);
f1m=nodalf(s,u-del*[qr+r;qz]); f2m=nodalf(s,u-del*[qr-r;qz]);
sig1=0.25*pr'*(f1m+f1p-f2m-f2p);
f1p=nodalf(s,u+del*[qi+r;qz]); f2p=nodalf(s,u+del*[qi-r;qz]);
f1m=nodalf(s,u-del*[qi+r;qz]); f2m=nodalf(s,u-del*[qi-r;qz]);
sig2=0.25*pr'*(f1m+f1p-f2m-f2p);
sig=sig1+sig2; % sig, pause
% step 5
f1p=nodalf(s,u+del*[qr+sr;qz]); f2p=nodalf(s,u+del*[qr-sr;qz]);
f1m=nodalf(s,u-del*[qr+sr;qz]); f2m=nodalf(s,u-del*[qr-sr;qz]);
d1=0.25*pr'*(f1m+f1p-f2m-f2p);
f1p=nodalf(s,u+del*[qi+si;qz]); f2p=nodalf(s,u+del*[qi-si;qz]);
f1m=nodalf(s,u-del*[qi+si;qz]); f2m=nodalf(s,u-del*[qi-si;qz]);
d2=0.25*pr'*(f1m+f1p-f2m-f2p);
f1p=nodalf(s,u+del*[qr+si;qz]); f2p=nodalf(s,u+del*[qr-si;qz]);
f1m=nodalf(s,u-del*[qr+si;qz]); f2m=nodalf(s,u-del*[qr-si;qz]);
d3=0.25*pr'*(f1m+f1p-f2m-f2p);
f1p=nodalf(s,u+del*[qi+sr;qz]); f2p=nodalf(s,u+del*[qi-sr;qz]);
f1m=nodalf(s,u-del*[qi+sr;qz]); f2m=nodalf(s,u-del*[qi-sr;qz]);
d4=0.25*pr'*(f1m+f1p-f2m-f2p);
d0=d1+d2+d3+d4; % d0, pause
% step 6, 3rd derivatives 
fdsw=2; % order of the FD 
switch fdsw; 
  case 2;  % seems sufficient, and consistent with other FDs
d=[qr;qz]; fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d); 
g1=pr'*(-0.5*fmm+fm-fp+0.5*fpp)/(del^3); 
d=[qi;qz]; fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d); 
g2=pi'*(-0.5*fmm+fm-fp+0.5*fpp)/(del^3); 
d=[qr+qi;qz]; fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d); 
g3=(pr+pi)'*(-0.5*fmm+fm-fp+0.5*fpp)/(del^3); 
d=[qr-qi;qz]; fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d); 
g4=(pr-pi)'*(-0.5*fmm+fm-fp+0.5*fpp)/(del^3); 
    case 4; % 
d=[qr;qz]; fmmm=nodalf(s,u-3*del*d); fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d);  fppp=nodalf(s,u+3*del*d); 
g1=pr'*(fmmm/8-fmm+13*fm/8-13*fp/8+fpp-fppp/8)/(del^3); 
d=[qi;qz]; fmmm=nodalf(s,u-3*del*d); fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d);  fppp=nodalf(s,u+3*del*d); 
g2=pi'*(fmmm/8-fmm+13*fm/8-13*fp/8+fpp-fppp/8)/(del^3); 
d=[qr+qi;qz]; fmmm=nodalf(s,u-3*del*d); fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d);  fppp=nodalf(s,u+3*del*d); 
g3=(pr+pi)'*(fmmm/8-fmm+13*fm/8-13*fp/8+fpp-fppp/8)/(del^3); 
d=[qr-qi;qz]; fmmm=nodalf(s,u-3*del*d); fmm=nodalf(s,u-2*del*d); fm=nodalf(s,u-del*d); 
fp=nodalf(s,u+del*d); fpp=nodalf(s,u+2*del*d);  fppp=nodalf(s,u+3*del*d); 
g4=(pr-pi)'*(fmmm/8-fmm+13*fm/8-13*fp/8+fpp-fppp/8)/(del^3); 
end
g0=2*(g1+g2)/3+(g3+g4)/6; % g0, pause
%fprintf('%g, %g, %g, %g, %g\n', g1,g2,g3,g4,g0); 
% step 7
ga=0.5*(g0-2*sig+d0)/abs(om); % l1(0) according to Kuz04 (p536) 
c1=abs(om)*ga; % Kuz04, (3.19) 
% compute real(mu'(lam_H))
p=s; lam=getlam(p); lamp=lam+dellam; 
p=setlam(p,lamp); u=nloop(p,p.u); r=resi(p,u); Gu=getGu(p,u,r); 
[ineg,muv,V]=spcalc(Gu,p,p.sol.j0); mu=real(muv(1)); 
mup=-(mu-mu0)/dellam; %mup, c1
if ((c1>0 && mup>0) || (c1<0 && mup<0)); dlam=-1;  % compute direction of lam 
else dlam=1; end
al=sqrt(-dlam*mup/c1); % amplitude coefficient 
end