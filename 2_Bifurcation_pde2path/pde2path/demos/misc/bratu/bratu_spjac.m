function Gvvph=bratu_spjac(p,u) 
% get partial_u (G_u phi), called in getGu if p.spcontsw==1 
par=u(2*p.nu+1:end); ph=u(p.nu+1:2*p.nu); u=u(1:p.nu); % params, Evec, PDE-vars
fuu=10*par(1)*exp(u); Gvvph=-p.mat.M*diag(fuu.*ph); 
end