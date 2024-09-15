function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; fem=p.pdeo.fem; n=p.np; ov=ones(1,n); 
po=getpte(p); th=po(2,:)'; dd=cos(th); 
%[Kphi,M,~]=fem.assema(gr,[ov 0*ov; 0*ov 0*ov],1,1); [Kth,~,~]=fem.assema(gr,[0*ov 0*ov;0*ov dd'],1,1); 
ov=ones(n,1); [Kphi,M,~]=fem.assema(gr,[ov; 0*ov],1,1); % old 
[Kth,~,~]=fem.assema(gr,[0*ov;dd],1,1); 
LB=spdiags(1./dd,0,n,n)*Kth+spdiags((1./dd).^2,0,n,n)*Kphi; 
p.mat.K=filltrafo(p,LB); p.mat.M=filltrafo(p,M); 
Dphi=convection(fem,gr,[1; 0]); p.mat.Dphi=filltrafo(p,Dphi); 