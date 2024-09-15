%% init and zero-branch 
p=[]; lx=sqrt(2)*pi/2; ly=lx; lz=8*lx; ndim=3; lam=-0.001; nu=1.5; par=[lam; nu];  
%nx=8; sw.sym=1; sw.ref=0;  dir='t1'; % np=17375, slow! 
nx=6; sw.sym=1; sw.ref=0;  dir='t1'; % np=6923, nice. 
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); figure(1); p.pdeo.grid.plotFaces; p.np 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; p=setfn(p,dir); 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; z=po(3,:)'; % extract coord from p 
xmin=min(x); ymin=min(y); x=x-xmin; y=y-ymin;  
p.sw.spcalc=0; p.sw.bifcheck=0;  p.plot.pstyle=3; 
%% BCC2zero from iguess and Newton; note that we only seed u1, while u2=0 
s2=sqrt(2); zc=0; lam=-0.3; nu=1.5; 
bca=0.3; ta=0; la=0; % bcc-amplitude, tube-ampl. and lamella-ampl. for guess
u=la*cos(z/s2).*(z<=zc); u=u+ta*((cos(x+y)/s2)+cos((x-y)/s2)).*(z<=zc);  
u=u+bca*((cos((x+y)/s2)+cos((x-y)/s2))+(cos((x+z)/s2)+cos((x-z)/s2))).*(z>zc); 
u=u+bca*(cos((y+z)/s2)+cos((y-z)/s2)).*(z>zc); v=[30 45]; 
p.u(1:p.np)=u; p.u(p.nu+1)=lam; p.u(p.nu+2)=nu;  r=norm(resi(p,p.u),'inf'); 
spl(p,''); view(v);  title('guess 1'); p=setfn(p,'b2z'); p.nc.imax=15; p.file.smod=1; 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu);
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); spl(p,''); 
huclean(p); p=resetc(p); p.sol.restart=1; p.nc.imax=5; p.sol.ds=0.01; p=pmcont(p,400);
% ptli=[1 10]; aux.cfile=0; inegli=fspcalc('b2z',ptli,aux)  % to a-posteriori compute stability 
%% BCC2tubes front from iguess and Newton 
s2=sqrt(2); zc=pi; bca=0.4; ta=1; la=0; lam=0.2; nu=1.5; % bcc2tubes, works with lz=8*lx
u=la*cos(z/s2).*(z<=zc); u=u+ta*(cos((x+y)/s2)+cos((x-y)/s2)).*(z<=zc);  
u=u+bca*((cos((x+y)/s2)+cos((x-y)/s2))+(cos((x+z)/s2)+cos((x-z)/s2))).*(z>zc); 
u=u+bca*(cos((y+z)/s2)+cos((y-z)/s2)).*(z>zc);
p.u(1:p.np)=u; p.u(p.nu+1)=lam; p.u(p.nu+2)=nu; r=norm(resi(p,p.u),'inf'); 
spl(p,''); view(v);  title('guess 2'); p=setfn(p,'b2t'); p.nc.imax=15; 
% t1=0; ts=[]; nc=0; dt=0.01; nt=1000; pmod=200; smod=500; p.mat.Kadv=0; % for time stepping
%% the tint loop, repeat this cell until residual is small (here just once) 
% [p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalf); % then plot time series of res 
%tss=2; figure(3); clf; plot(ts(1,tss:end),ts(2,tss:end)); axis tight; legend('res'); xlabel('t');
%% Newton
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); spl(p,''); 
%% cont of obtained solution
clf(2); p=resetc(p); p.sol.restart=1; p.nc.imax=5; p.sol.ds=0.01; p.file.smod=50; p=pmcont(p,400); 
%% generate (pure) BCC and tube branches over minimal domain  
p=[]; lx=sqrt(2)*pi/2; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=1.5; par=[lam; nu];  
nx=3; sw.sym=1; sw.ref=2; dir='BCC1s'; p.plot.pstyle=2; 
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p.pdeo.grid.plotFaces; p.np, pause 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; p=setfn(p,dir); p.plot.pstyle=3; 
p.sw.bifcheck=2; p.pm.resfac=1e-3; p.usrlam=[0 0.2]; p.nc.lammax=0.4; p=cont(p,10); 
%% qswibra; 3 tubes as kernel; afterwards set some switches 
aux=[]; aux.isotol=1e-16; aux.m=3; p0=qswibra('BCC1s','bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-4; p0.sw.spcalc=1; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.sw.foldcheck=0; p0.sw.spcalc=1; p0.plot.pstyle=2; 
%% select BCC and cont in hot direction, use gentau for tubes;
%p=seltau(p0,1,'sBCC',2); p.sol.ds=-0.02; p=pmcont(p,40); 
p=gentau(p0,1,'tu1'); p.sol.ds=-0.02; p=pmcont(p,30); 
%% plot BD 
fnr=3; cmp=3; figure(fnr); clf; 
plotbra('sBCC','pt33',fnr,cmp,'cl','r','lab',29); 
plotbra('tu1','pt30',fnr,cmp,'cl','b','lab',26); 
plotbra('b2z','pt360',fnr,cmp,'cl',p2pc('r1'),'lp',350,'lab',[0, 200, 350],'lw',1); 
plotbra('b2t','pt350',fnr,cmp,'cl',p2pc('r2'),'lab',[0 100 350],'lw',1); 
xlabel('\lambda'), ylabel('||u||'); box on; 
%% soln plots 
v=[-60 45]; 
plotsol('sBCC','pt29',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); colorbar 'off'; pause
plotsol('tu1','pt26',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); colorbar 'off'; pause
v=[30 45]; plotsol('b2z','pt0',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); pause 
plotsol('b2z','pt350',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); 
%% 
ps=2; for i=[0 100 350]
plotsol('b2t',['pt' mat2str(i)],1,1,ps); view(v); xlabel(''); ylabel(''); zlabel(''); axis image; pause
end 