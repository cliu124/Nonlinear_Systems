function out=secobraHPC(p,u)  % mod of hobra.m to compute deviation from tr. state 
switch p.sw.spcont
    case 0; np=p.nu/p.nc.neq; n=p.nu; % use nu, i.e. all compos of u 
    case 1; np=p.nu/p.nc.neq; n=p.nu/2;
    case 2; np=p.nu/p.nc.neq; n=p.nu/2;
    case 3; np=p.nu/p.nc.neq; n=p.nu/3;
end
par=u(p.nu+1:end); [as,us]=getss(par); 
if p.sw.para>2 % Hopf setup; get data from p.hopf
    ho=p.hopf; m1=ho.T; 
    m2=max(max(ho.y(1:np,:)-as)); m3=min(min(ho.y(1:np,:)-as)); 
    l2v=zeros(1,ho.tl); 
    for i=1:length(ho.t); l2v(i)=(ho.y(1:np,i)'-as)*(ho.tom.M(1:np,1:np)*(ho.y(1:np,i)-as)); end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/p.vol); 
else % steady state continuation  
    m1=0; m2=max(u(1:np)-as); m3=min(u(1:np)-as); 
    l2=(u(1:np)-as)'*(p.mat.M(1:np,1:np)*(u(1:np)-as)); m4=sqrt(l2/p.vol); 
end
%p.nu, p.nc.nq
try; om=u(p.nu+p.nc.nq+5); catch; om=0; end 
out=[par; m1; m2; m3; m4; om]; 