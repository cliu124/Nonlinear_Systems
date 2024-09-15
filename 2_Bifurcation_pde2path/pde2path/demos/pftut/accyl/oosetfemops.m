function p=oosetfemops(p) % AC on cylinder with lid; two fields u1,u2
gr=p.pdeo.grid; fem=p.pdeo.fem; n=p.np; 
[Kphi,M,~]=fem.assema(gr,[ones(n,1);0*ones(n,1)],1,1); 
[Kz,~,~]=fem.assema(gr,[0*ones(n,1); ones(n,1)],0,1); 
p.nu=p.nu1; % set nu to nu1 (dim of first field) for filltrafo 
p.mat.Kphi=filltrafo(p,Kphi); p.mat.Kz=filltrafo(p,Kz); 
p.mat.M0=M; p.mat.M1=filltrafo(p,M); 
Dphi=convection(fem,gr,[1;0]); p.mat.Dphi=filltrafo(p,Dphi); % for PC 
p.nu=p.nu1+p.np2; % reset nu to total number of unknowns
gr=p.p2.grid; fem=p.p2.fem; % assemble matrices for 2nd field 
[K,M,~]=fem.assema(gr,1,1,1); p.mat.K2=K; p.mat.M2=M; 
% put together full mass matrix M needed for spcalc 
p.mat.M=[p.mat.M1 sparse(p.nu1,p.np2); sparse(p.np2,p.nu1) M]; 