function out=gpbra3(p,u) % GPbra 
nu=p.nu; np=p.np; par=u(nu+1:end); u1=u(1:nu); 
N=sum(p.mat.M*(u1.^2)); d=soldiff(p,u); %H=ham(p,u); 
out=[par; max(abs(u1(1:np))); N; d]; 
end

function d=soldiff(p,u)
par=u(p.nu+1:end); u=u(1:p.nu); lam=par(1); s=par(2); sig=par(3); ga=par(4); 
x=getpte(p); x=x'; bet=(-ga/lam)^(1/(2*sig)); V=pot(x,s); 
unls=(1/bet)*(sig+1)^(1/(2*sig))./(cosh(sqrt(lam)*sig*(x-s)).^(1/sig));  
figure(1); hold on; plot(x,unls,x,V); hold off;
d=sqrt(sum(p.mat.M*((u-unls).^2)));
end 
