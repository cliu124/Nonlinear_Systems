function out=hcylbra(p,u)
par=u(p.nu+1:end); l1=par(2); c0=par(3); np=p.np; 
Ai=getA(p,u); V=getV(p,u); M=getM(p); 
M=M(1:np,1:np); H=p.u(np+1:2*np); M(p.idx,p.idx)=0; % ignore bdries for E 
E=sum(M*(H-c0).^2)+l1*Ai; % energy 
mqd=meshqdat(p); % max(A),min(A), max(h/r), max(h/R), max(h), min(h) 
out=[u(p.nu+1:end); Ai; V; E; l1*Ai; mqd]; 