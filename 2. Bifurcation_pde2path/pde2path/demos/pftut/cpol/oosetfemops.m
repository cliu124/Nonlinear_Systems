function p=oosetfemops(p) % AC on sphere with coupling to bulk 
gr=p.pdeo.grid; fem=p.pdeo.fem; n=p.nps; 
po=getpte(p); th=po(2,:)'; dd=cos(th); ov=ones(n,1); 
[Kphi,M,~]=fem.assema(gr,[ov; 0*ov],dd,1); [Kth,~,~]=fem.assema(gr,[0*ov;dd],1,1); 
LB=Kth+spdiags(1./dd,0,n,n)*Kphi; 
p.nu=p.nus; % set nu to nus (dim of first field) for filltrafo 
p.mat.K=filltrafo(p,LB); p.mat.M1=filltrafo(p,M); 
Dphi=convection(fem,gr,[1;0]); p.mat.Dphi=filltrafo(p,Dphi); % for PC 
p.nu=p.nus+p.npb; % reset nu to total number of unknowns
gr=p.p2.grid; fem=p.p2.fem; % assemble matrices for bulk field 
[K,M,~]=fem.assema(gr,1,1,1); p.mat.Kb=K; p.mat.Mb=M; 
% assemble full mass M needed for spcalc 
p.mat.M=[p.mat.M1 sparse(p.nus,p.nub); sparse(p.nub,p.nus) p.mat.Mb]; 
bc=gr.robinBC(1,1); % BCs for the 2nd set of BC matrices
gr.makeBoundaryMatrix(bc); 
[p.mat.Q2,p.mat.G2,~,~]=fem.assemb(gr);
p.mat.vM1=sum(p.mat.M1,1); % \int_\Ga u1 dS=R^2*vM1*u1
p.mat.vM2=sum(p.mat.Mb,1); % \int_\Om u2 dx=vM2*u2