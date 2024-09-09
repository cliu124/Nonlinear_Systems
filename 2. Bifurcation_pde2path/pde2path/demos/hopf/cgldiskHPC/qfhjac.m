function q_u=qfhjac(p,y) % u_derivatives of qfh
tl=p.hopf.tl; n=p.nu; 
q_u=zeros(1,tl*p.nu);
u0x=p.u0x(1:p.nu);
for i=1:tl; q_u((i-1)*p.nu+1:i*p.nu)=u0x'; end % Jac of phase constr.