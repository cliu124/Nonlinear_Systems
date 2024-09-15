%%
% hobra: output for bifurcation diagram, setup for Hopf
% [T,max(u),min(u),||u||]
function out=hobra(p,u)  
gotit=0; np=p.nu/p.nc.neq; M=p.mat.M(1:np,1:np);
if isfield(p,'hopf'); ho=p.hopf; 
    if isfield(ho,'y'); m1=ho.T; 
      m2=max(max(ho.y(1:np,:))); m3=min(min(ho.y(1:np,:))); 
      l2v=zeros(1,ho.tl); 
      for i=1:length(ho.t); l2v(i)=ho.y(1:np,i)'*(ho.tom.M(1:np,1:np)*ho.y(1:np,i)); end 
      l2=trapz(ho.t,l2v); m4=sqrt(l2/(p.vol)); 
      m5=sum(M*(ho.y(1:np,1)+ho.y(np+1:2*np,1))); 
      m6=max(max(ho.y(np+1:2*np,:))); m7=min(min(ho.y(np+1:2*np,:))); 
    gotit=1; end; end 
if ~gotit; m1=0; m2=max(u(1:np)); m3=min(u(1:np)); 
    l2=u(1:np)'*(p.mat.M(1:np,1:np)*u(1:np)); m4=sqrt(l2/p.vol); 
    m5=sum(M*(u(1:np)+u(np+1:2*np))); m6=max(u(np+1:2*np)); m7=min(u(np+1:2*np));
end
out=[u(p.nu+1:end); % parameters 
        m1; m2; m3; m4; m5; m6; m7]; 
    