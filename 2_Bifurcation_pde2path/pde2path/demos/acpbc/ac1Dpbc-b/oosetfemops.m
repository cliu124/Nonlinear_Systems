function p=oosetfemops(p) % for 1Dpbc, O(2) equiv, hence need phase-condition 
% <u,\pa_x u0>=0, where u0=reference profile. 
% Thus need Kx with M\pa_x u0=Kx*u0
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[K,M,~]=fem.assema(gr,1,1,1); 
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
Kx=convection(fem,gr,1); p.mat.Kx=filltrafo(p,Kx); 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % transform of K and M 