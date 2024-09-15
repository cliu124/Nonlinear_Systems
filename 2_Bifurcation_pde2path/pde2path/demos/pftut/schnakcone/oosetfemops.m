function p=oosetfemops(p) % Schnak on cone 
par=p.u(p.nu+1:end); a=par(4); [K,M]=LBcone(p,a); 
M=[M 0*M; 0*M M]; p.mat.M=M; p.mat.K=K;
