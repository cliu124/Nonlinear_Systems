function mu=ploteigsys(dir,fname,m)
% ploteigsys: convenience function to plot m evecs, starting with most negative 
if ischar(dir); p=loadp(dir,fname);
else; p=dir; if nargin>1; m=fname; end % input is struct, m is 2nd argument 
end 
%p.nc.eigref=-1000; 
fss=p.plot.fs; p.plot.fs=12; u=p.u; lam=getlam(p);
r=resi(p,u); [Gu,Glam]=getder(p,u,r);  % residual and jacobian
M=getM(p); %p.mat.M; 
if(p.nc.nq>0) % trivially extend M in case of auxiliary equations
    if isfield(p.fuha,'qMfu'); [qL,qU,qD]=p.fuha.qMfu(p); 
    else [qL,qU,qD]=stanqM(p); 
    end 
    M=getM(p); M=[[M(1:p.nu,1:p.nu) qU]; [qL qD]]; 
end
vs=size(Gu,1); opts.v0=ones(vs,1)/vs; opts.disp=0; pstyle=p.plot.pstyle; 
[phiv,mu]=myeigs(Gu,M,m,p.nc.eigref,opts,p); 
phiv=real(phiv);   % getting EVecs
mud=diag(mu); [muds,idx]=sort(real(mud)); mu=mud(idx)';
phiv=phiv(:,idx); 
fprintf(['lam=' num2str(lam), '; ' num2str(m) ' smallest eigenvalues: ' num2str(mu(1:m),2) '\n']);
try Xcont=p.sw.Xcont; catch Xcont=0; end 
par=p.u(p.nu+1:end); 
set(0,'Units','pixels'); scnsize = get(0,'ScreenSize'); 
scnw=scnsize(3); scnh=scnsize(4); sz=scnh/3-scnh*5/100; 
for i=1:m
   switch mod(i,3); case 1; pheight=2*scnh/3; case 2; pheight=sz+115; case 0; pheight=0; end;
   pw=ceil(i/3); pw=pw+1; set(figure(10+i),'Position',[scnw-pw*sz pheight sz sz]);
   if Xcont>0; p.up=[phiv(1:p.nu,i);p.u(p.nu+1:end)]; end 
   plotsolu(p,[phiv(:,i);par],10+i,1,pstyle); 
   title([p.file.pname mat2str(p.file.count-1) ', \mu_1=' mat2str(mu(i),4) ', \phi_' mat2str(i,2) ':']); 
end 
drawnow; 