function p=oosetfemops(p) % for brusselator with time periodic forcing  
[p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.Ms=M; p.mat.M0=[M 0*M; 0*M M]; % PDE-part M 
p.mat.M=[[p.mat.M0 zeros(2*p.np,2)]; % extend PDE-M by oscillator M
         [zeros(2,2*p.np) speye(2)]]; 