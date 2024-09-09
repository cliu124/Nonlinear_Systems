function r=pderesi(p,u)
% PDERESI(p,u): return pde-part of residual
if p.sw.sfem==0 % residual with full assembling
    [c,a,f,b]=p.fuha.G(p,u); bc=p.fuha.bc(p,u); upde=p.mat.fill*u(1:p.nu); 
    [K,F]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,c,a,f,upde); 
    if(any(b)); Kadv=assemadv(p.mesh.p,p.mesh.t,b); K=K-Kadv; end 
    K=filltrafo(p,K); 
    F=p.mat.fill'*F; % adapt in case of periodic domain
    r=K*u(1:p.nu)-F;
else r=p.fuha.sG(p,u); % use "simple" assembling of G
end
