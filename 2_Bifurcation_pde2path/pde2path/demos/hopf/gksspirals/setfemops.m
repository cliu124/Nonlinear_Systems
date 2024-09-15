function p=setfemops(p) % mod of standard setfemops to compute the 
% rotation matrix Krot
if p.sw.sfem<0; p=oosetfemops(p); return; end
try upde=p.mat.fill*p.u(1:p.nu); % set to full domain vector
catch; upde=zeros(p.nu,1); end 
% fold/branch point continuation: set neq to original value:
if (p.sw.spcont~=0) neq=p.nc.neq/2; upde=upde(1:p.np*neq); else neq=p.nc.neq; end
bc=p.fuha.bc(p,p.u); 
[~,p.mat.M,~,~,p.mat.bcG,~,~]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,...
    0,1,zeros(neq,1),upde); p.mat.M=filltrafo(p,p.mat.M); 
if(p.sw.sfem==1) % also K and Kadv needed 
    [p.mat.K,~]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,p.eqn.c,p.eqn.a,...
      zeros(neq,1)); p.mat.K=filltrafo(p,p.mat.K);
    if any(p.eqn.b)
        p.mat.Kadv=assemadv(p.mesh.p,p.mesh.t,p.eqn.b); 
        p.mat.Kadv=filltrafo(p,p.mat.Kadv);
    else
        p.mat.Kadv=0; p.mat.Kadv=sparse(p.mat.Kadv); 
    end
end
po=getpte(p); x=po(1,:)'; y=po(2,:)'; % compute rot-matrix
xi=pdeintrp(p.mesh.p,p.mesh.t,x); yi=pdeintrp(p.mesh.p,p.mesh.t,y);
b=zeros(p.nc.neq^2,p.mesh.nt); b(1,:)=yi; b(2,:)=-xi;
Krot=-assemadv1(p.mesh.p,p.mesh.t,b); 
p.mat.Krot=[Krot 0*Krot; 0*Krot Krot]; 


