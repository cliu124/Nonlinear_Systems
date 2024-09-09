function Kvud=getKuvd(p,par,u,v)
% getKuvd: compute \pa_u (div(c(u)\grad v) via numjac; useful for quasilin problems 
% 1 component version, i.e., size(u)=size(v)=p.np
% see, e.g. demo /demos/acsuite/acql 
global pj vj parj; pj=p; vj=v; parj=par; 
n=length(u); thresh=p.nc.del*ones(n,1); pj=p; pj.u=u; 
M=getM(p); M=M(1:n,1:n); 
gr=p.pdeo.grid; fem=p.pdeo.fem;  ut=p.mat.p2c*u; 
cc=p.fuha.cfu(ut,par); [K,~,~]=fem.assema(gr,cc,0,0);  r=K*v; 
S=M>0; 
[Kvud,njfac,njG]=numjac('Kuv',0,u,r,thresh,[],0,S,[]); 