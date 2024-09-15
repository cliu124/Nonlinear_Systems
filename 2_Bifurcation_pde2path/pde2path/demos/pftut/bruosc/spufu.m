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
n=p.np;nu=p.nu;par=p.u(p.nu+1:end); u=[p.u(1); p.u(n+1)]; % short Hom-state vector 
u=[u;par]; p.np=1; p.nu=2; [f1u,f1v,f2u,f2v]=njac(p,u); J=[[f1u f1v]; [f2u f2v]]; 
kv=0:0.01:1; kl=length(kv); muv=zeros(2,kl); % provide wave-nrs and mem 
du=par(3); dv=par(4); % diffusion param.; this and k-range usually only problem dep. things 
for i=1:kl %  now loop over wave-nr and compute Evals 
    k=kv(i); K=[[du*k^2 0];[0 dv*k^2]];  A=J-K; % Jac in Fourier space 
    mu=eig(A); [mus, ix]=sort(real(mu),'descend'); % sorted eigenvalues 
    for j=1:2; muv(j,i)=mu(ix(j)); end 
end 
figure(10); clf; plot(kv, real(muv(1,:)), kv, imag(muv(1,:))); % plot leading Eval   
title(['par=' mat2str(par(1:4),3)],'fontsize',p.plot.fs); set(gca,'FontSize',p.plot.fs); 
legend('real','imag'); pause; p.np=n; p.nu=nu; % reset # spatial points to full problem 
end
