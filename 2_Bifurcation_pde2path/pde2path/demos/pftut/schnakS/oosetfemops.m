function p=oosetfemops(p)
gr=p.pdeo.grid; fem=p.pdeo.fem; par=p.u(p.nu+1:end); R=par(4); 
po=getpte(p); th=po(2,:)'; dd=cos(th); n=p.np; ov=ones(n,1); 
[Kphi,M,~]=fem.assema(gr,[ov; 0*ov],1,1); [Kth,~,~]=fem.assema(gr,[0*ov;dd],1,1); 
LB=spdiags(1./dd,0,n,n)*Kth+spdiags((1./dd).^2,0,n,n)*Kphi; 
p.mat.K=filltrafo(p,LB); 
Dphi=convection(fem,gr,[1;0]); Dphi=[Dphi 0*Dphi; 0*Dphi Dphi]; p.mat.Dphi=filltrafo(p,Dphi); 
M=[M 0*M; 0*M M]; p.mat.M=filltrafo(p,M);
