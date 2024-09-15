function q_u=qfhjac(p,y) % u_derivatives of qfh
M=p.mat.M0/p.vol; tl=p.hopf.tl; n=p.nu; 
q_u=zeros(2,tl*p.nu);
for i=1:1; q_u(1,(i-1)*p.nu+1:i*p.nu)=M*ones(n,1); end % Jac of mass constr.
if isfield(p,'u0x'); u0x=p.u0x(1:p.nu); else u0x=p.mat.Kx*p.u(1:p.nu); end 
for i=1:tl; q_u(2,(i-1)*p.nu+1:i*p.nu)=u0x'; end % Jac of phase constr.