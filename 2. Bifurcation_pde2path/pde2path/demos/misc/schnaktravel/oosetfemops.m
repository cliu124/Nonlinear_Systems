function p=oosetfemops(p) 
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); 
Kx=convection(p.pdeo.fem,gr,1); % x-derivative 
N=0*sparse(p.np,p.np); M=[[M N];[N M]]; p.mat.M=M; % system M 
d=60; p.mat.K=[[K N]; [N d*K]]; p.mat.Kx=[[Kx N]; [N Kx]];
if p.sw.bcper>0; [p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
 p.mat.K=filltrafo(p,p.mat.K); 
 p.mat.M=filltrafo(p,p.mat.M);
 p.mat.Kx=filltrafo(p,p.mat.Kx);  
end 


