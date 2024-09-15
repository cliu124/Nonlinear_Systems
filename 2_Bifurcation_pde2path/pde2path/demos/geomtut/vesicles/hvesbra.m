function out=hvesbra(p,u) % for closed vesicles
par=u(p.nu+1:end); l1=par(2); l2=par(3); c0=par(4); np=p.np;  
A=getA(p,p.u); V=getV(p,p.u); N=getN(p,p.X); M=getM(p,p.X+u(1:p.np).*N); 
M=M(1:np,1:np); H=u(np+1:2*np); M(p.idx,p.idx)=0; 
E0=sum(M*(H-c0).^2)/(16*pi); % bending energy 
E1=E0-l1/2*A+l2*V;  % Lagr.
mqd=meshqdat(p); % q; max(A); min(A); max(h); min(h)
A=par(5);
R=sqrt(A/(4*pi)); V0=4*pi*R/3; v=V/V0; c=c0*R; V1=altvol(p); 
Gb=sum(M*(H.^2))/(16*pi); Hi=sum(M*H); del_a=Hi/(4*pi*R); % data for BiC model
out=[par;  A;  V; E0; E1; l1*A; v;  c; V1; mqd; Gb; Hi; del_a]; 
%    1     13         16        18         21-25 26  27   28
% mqd(1)=del used in refufu.m for mesh refinement. 
% del_a can be used to plot branch in BiC-Model sense 