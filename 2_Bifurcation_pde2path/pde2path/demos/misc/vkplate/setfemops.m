function p=setfemops(p) % modified for vkplate: generate additional K2 which 
% will be multiplied by lam in vksG and vksGjac 
% (Note that for us with fold/branch point continuation need to adjust vecto
% dimensions as in setfemops from teh p2plib.)
upde=p.mat.fill*p.u(1:p.nu); neq=p.nc.neq; m=p.mesh;
bc=p.fuha.bc(p,p.u); 
[~,p.mat.M,~,~,p.mat.bcG,~,~]=assempde(bc,m.p,m.e,m.t,0,1,zeros(neq,1),upde);
[p.mat.K,~]=assempde(bc,m.p,m.e,m.t,p.eqn.c,p.eqn.a, zeros(neq,1),upde);
[p.mat.K2,~]=assema(p.mesh.p,p.mesh.t,p.eqn.c2,0,zeros(neq,1));


