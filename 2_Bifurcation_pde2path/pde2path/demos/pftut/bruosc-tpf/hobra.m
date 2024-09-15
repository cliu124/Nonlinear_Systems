function out=hobra(p,u)  
% hobra: output for bifurcation diagram, setup for Hopf
% [paras, T,max(u),min(u),||u||]
switch p.sw.spcont
    case 0; np=p.nu/p.nc.neq; n=p.nu; % use nu, i.e. all compos of u 
    case 1; np=p.nu/p.nc.neq/2; n=p.nu/2;
    case 2; np=p.nu/p.nc.neq/2; n=p.nu/2;
    case 3; np=p.nu/p.nc.neq/3; n=p.nu/3;    
end
par=u(p.nu+1:end); % par', pause 
if p.sw.para>2 % Hopf setup; get data from p.hopf
    ho=p.hopf; m1=ho.T; m2=max(max(ho.y(1:n,:))); m3=min(min(ho.y(1:n,:))); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=ho.y(1:n,i)'*(ho.tom.M(1:n,1:n)*ho.y(1:n,i)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
    m5=max(max(ho.y(p.nu-1,:))); m6=min(min(ho.y(p.nu,:))); 
else % steady state continuation  
    m1=u(p.nu); m2=max(u(1:n)); m3=min(u(1:n)); 
    l2=u(1:n)'*(p.mat.M(1:n,1:n)*u(1:n)); m4=sqrt(l2/p.vol); 
    m5=0; m6=0; 
end
out=[par; m1; m2; m3; m4; m5; m6]; 
    