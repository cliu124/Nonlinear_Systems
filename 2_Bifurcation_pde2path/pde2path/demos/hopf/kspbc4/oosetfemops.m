function p=oosetfemops(p) % with filltrafo to transform to per.domain
gr=p.pdeo.grid; [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); Kx=convection(p.pdeo.fem,gr,1);
p.mat.K=filltrafo(p,K); M=filltrafo(p,M); p.mat.Kx=filltrafo(p,Kx);
p.mat.M0=M; p.mat.M=M^2; % save M as M0 and redefine M for the 4th order setup