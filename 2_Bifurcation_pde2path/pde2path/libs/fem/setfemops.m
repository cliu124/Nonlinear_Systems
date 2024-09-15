function p=setfemops(p)
% SETFEMOPS: generate and store FEM operators 
% Uses full domain assembling, and transforms in case of periodic bc.
%
%  p=setfemops(p)
%
% See also loadp, rec2cyl, rec2tor, filltrafo, assemadv, assempde
if p.sw.sfem<0; 
   if isfield(p.fuha,'setops'); p=p.fuha.setops(p); % possible function handle 
   else p=oosetfemops(p); % for backward compatibility 
   end
   return; 
end
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

