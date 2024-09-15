function out=hobra(p,u)  
% hobra: output for bifurcation diagram, setup for Hopf, including Krylov 
% [paras, T,max(u),min(u),||u||]
switch p.sw.spcont
    case 0; np=p.nu/p.nc.neq; % only use 1st compo (overload if wished differently!) 
    case 1; np=p.nu/p.nc.neq; 
    case 2; np=p.nu/p.nc.neq; 
    case 3; np=p.nu/p.nc.neq; 
end
par=u(p.nu+1:end); 
if p.sw.para<3 %steady state continuation  
    m1=0; m2=max(u(1:np)); m3=min(u(1:np)); 
    l2=u(1:np)'*(p.mat.M(1:np,1:np)*u(1:np)); m4=sqrt(l2/p.vol); 
else % Hopf setup; 
    if p.sw.para<6;   % legacy, with t-dependence 
    ho=p.hopf; m1=ho.T; 
    m2=max(max(ho.y(1:np,:))); m3=min(min(ho.y(1:np,:))); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=ho.y(1:np,i)'*(ho.tom.M(1:np,1:np)*ho.y(1:np,i)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
    else % shooting  
    m1=u(p.nu+p.hopf.iT); par=p.u(p.nu+1:end); 
    [uv,ts]=tintup(p,m1,p.hopf.nt,u(1:p.nu),par); ts=ts/m1;  %ts, pause 
    m2=max(max(uv(1:p.np,:))); m3=min(min(uv(1:p.np,:))); 
    tl=length(ts); l2v=0*ts;  % ts, np, size(ts), size(uv), tl, pause 
    for i=1:tl; l2v(i)=uv(1:np,i)'*(p.mat.M(1:np,1:np)*uv(1:np,i)); end 
    l2=trapz(ts,l2v); m4=sqrt(l2/p.vol); 
    %l2=u(1:np)'*(p.mat.M(1:np,1:np)*u(1:np)); m4=sqrt(l2/p.vol);     
    end    
end
out=[par; m1; m2; m3; m4]; 
    