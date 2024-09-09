function p=oosetfemops(p) % for 1Dpbc, with x-dep. K 
x=getpte(p); c=1+0.1*x.^2; [K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,c,1,1); 
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % transform of K and M 