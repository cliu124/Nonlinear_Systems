function r=pderesi(p,u)
% compute pde-part of residual
if p.sw.sfem~=1 % residual with full assembling
    upde=p.mat.fill*p.u(1:p.nu); 
    [c,a,f,b]=p.fuha.G(p,u); bc=p.fuha.bc(p,u); 
    [K,F]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde); 
    if(any(b)); Kadv=assemadv(p.mesh.p,p.mesh.t,b); K=K-Kadv; end 
    K=filltrafo(p,K); 
    F=p.mat.fill'*F; % adapt in case of periodic domain
    r=K*u(1:p.nu)-F;
else r=p.fuha.sG(p,u); % use "simple f"
end
if p.pffac>0; r(p.nu/2+p.pfn)=p.pffac*u(p.nu/2+p.pfn); % u2(pn)=0 
end


