function Kuv=Kuv(t,u) 
% Kuv: compute div(c(u) grad v) with c in p.fuha.cfu 
%
% used in getKuvd for \pa_u (div(c(u)\grad v) via numjac
global pj vj parj; 
p=pj; par=parj; n=length(u); gr=p.pdeo.grid; fem=p.pdeo.fem;  
ut=p.mat.p2c*u; cc=p.fuha.cfu(ut,par); 
[K,~,~]=fem.assema(gr,cc,0,0); 
Kuv=K*vj; 
end 