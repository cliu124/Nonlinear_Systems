function [c,a,f,b]=chemG(p,u)
% pde for chemotaxis model 
% u_t=0.25*Lap u-lam*div(u*grad v)+r*u(1-u)
% v_t=Lap v+(u/(1+u)-v) 
lam=u(p.nu+1:end); u=u(1:p.nu);
u=pdeintrp(p.mesh.p,p.mesh.t,u);duv=ones(1,p.mesh.nt); D=0.25; r=1.52;
f1=r*u(1,:).*(1-u(1,:)); f2=u(1,:)./(1+u(1,:))-u(2,:); a=0; b=0;
f=[f1;f2];
%c=isoc([[D*duv -lam*u(1,:)]; [0*duv duv]],p.neq,p.nt); return 
% set c by hand:
c1111=D*duv; c1112=0*duv; c1121=0*duv; c1122=D*duv; 
c1211=-lam*u(1,:); c1212=0*duv; c1221=0*duv;c1222=-lam*u(1,:);
c2111=0*duv; c2112=0*duv; c2121=0*duv; c2122=0*duv; 
c2211=duv; c2212=0*duv; c2221=0*duv; c2222=1*duv; 
c=[c1111; c1121; c1112; c1122; c2111; c2121; c2112; c2122; ...
   c1211; c1221; c1212; c1222; c2211; c2221; c2212; c2222];  