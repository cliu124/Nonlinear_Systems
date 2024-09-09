function speedplot(p,varargin)
tl=p.hopf.tl; nu=p.nu; sv=zeros(1,tl); 
par=p.u(p.nu+1:end); y=p.hopf.y; s=par(p.hopf.ilam); 
for i=1:tl; u=[y(1:p.nu,i);par];
    ux=p.mat.Kx*u(1:nu); utx=p.u0x; n=utx'*ux; 
    r=pderesi(p,u); G=r+s*ux; 
    sj=utx'*G/n; sv(i)=sj; 
end
lt='b'; if nargin>1; lt=varargin{1}; end 
figure(10); plot(p.hopf.T*p.hopf.t,sv,lt);  axis tight; 
set(gca,'FontSize',p.plot.fs); xlabel('t'); ylabel('s'); 