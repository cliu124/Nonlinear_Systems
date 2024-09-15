function out=hobratw(p,u)  
% hobratw: output for bifurcation diagram, setup for TW-Hopf
% [paras, T,max(u),min(u),||u||]
switch p.sw.spcont
    case 0; n=p.nu/p.nc.neq; % n=p.nu; % use nu, i.e. all compos of u 
    case 1; n=p.nu/p.nc.neq; %n=p.nu/2;    % BP-cont 
    case 2; n=p.nu/p.nc.neq; %n=p.nu/2;    % FP-cont 
    case 3; n=p.nu/p.nc.neq; %n=p.nu/3;    % HP-cont 
end
par=u(p.nu+1:end); % par', pause 
if p.sw.para>2 % Hopf setup; get data from p.hopf
    ho=p.hopf; m1=ho.T; m2=max(max(ho.y(1:n,:))); m3=min(min(ho.y(1:n,:))); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=ho.y(1:n,i)'*(ho.tom.M(1:n,1:n)*ho.y(1:n,i)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
else % steady state continuation  
    po=getpte(p); try L=p.L; catch; L=abs(po(1,1)-po(1,end)); end% domain size 
    m1=L/(p.kwnr*abs(par(p.spar))); % period for TW-cont (1D) 
    m2=max(u(1:n)); m3=min(u(1:n)); 
    l2=u(1:n)'*(p.mat.M(1:n,1:n)*u(1:n)); m4=sqrt(l2/p.vol); 
end
out=[par; m1; m2; m3; m4]; 
    