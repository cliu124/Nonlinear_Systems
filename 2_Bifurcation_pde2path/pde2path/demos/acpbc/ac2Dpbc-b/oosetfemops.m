function p=oosetfemops(p) % in problem-dir 
fem=p.pdeo.fem; gr=p.pdeo.grid; [K,M,~]=fem.assema(gr,1,1,1); % indep. of BC 
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % standard transforms of 
Kx=convection(fem,gr,[1;0]); Ky=convection(fem,gr,[0;1]); 
p.mat.Kx=filltrafo(p,Kx); p.mat.Ky=filltrafo(p,Ky); 