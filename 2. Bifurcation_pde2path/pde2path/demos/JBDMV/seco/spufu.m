%%
% spufu, adadption of STANUFU to plot spectrum of lin. around hom sol 
function [p,cstop]=spufu(p,brout,ds)
p=ulamcheck(p); % test if a desired lambda value has been passed
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
    p.file.count,printcon(getlam(p)),printcon(brplot),p.sol.res,p.sol.iter, ...
    p.sol.meth,printcon(ds));
if(p.sw.errcheck>0) fprintf('%5.2e ',p.sol.err);end
if(p.sw.spcalc==1) fprintf(' %2i ', p.sol.ineg); end;
npr=length(p.sw.bprint);
for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(p))+p.sw.bprint(i)))); end;
fprintf('\n');
cstop=0;
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end
% addition to STANUFU: plot the dispersion relation! 
n=p.np;nu=p.nu; par=p.u(p.nu+1:end);u=[p.u(1);p.u(n+1)]; % hom-state 2-vector 
u=[u;par];p.np=1;p.nu=2;[f1u,f1v,f2u,f2v]=nodaljac(p,u); J=[[f1u f1v];[f2u f2v]]; 
d=par(4); % diffusion param.; this and k-range usually only problem dep. things 
kmax=0.2; kv=0:0.002:kmax; kl=length(kv); muv=zeros(2,kl); % provide wave-nrs and mem 
x=getpte(p); lx=max(x); dk=pi/2/lx; 
kv2=0:dk:kmax; kl2=length(kv2); muv2=zeros(2,kl2);
for i=1:kl %  now loop over wave-nr and compute Evals 
    k=kv(i); K=[[k^2 0];[0 d*k^2]];  A=J-K; % Jac in Fourier space 
    mu=eig(A); [mus, ix]=sort(real(mu),'descend'); % sorted eigenvalues 
    for j=1:2; muv(j,i)=mu(ix(j)); end 
end 
for i=1:kl2 %  now loop over wave-nr and compute Evals 
    k=kv2(i); K=[[k^2 0];[0 d*k^2]];  A=J-K; % Jac in Fourier space 
    mu=eig(A); [mus, ix]=sort(real(mu),'descend'); % sorted eigenvalues 
    for j=1:2; muv2(j,i)=mu(ix(j)); end 
end 
figure(10); clf; plot(kv,real(muv(1,:)),kv,imag(muv(1,:)),kv,0*kv,'--k'); % plot leading Eval   
title(['par=' mat2str(par,3)],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
 p.np=n; p.nu=nu; % reset # spatial points to full problem 
hold on
plot(kv2,real(muv2(1,:)),'*b',kv2,imag(muv2(1,:)),'*r'); 
legend('real','imag'); xlabel('k'); title(''); axis tight; 
