% SH on long and slender bar, with localized BCCs from tint (DNS)  
% init and zero-branch 
p=[]; lx=sqrt(2)*pi/2; ly=lx; lz=8*lx; ndim=3; lam=-0.001; nu=1.5; par=[lam; nu];  
%nx=8; sw.sym=1; sw.ref=0;  dir='BCClong'; % np=17375, slow! 
nx=6; sw.sym=1; sw.ref=0;  dir='BCClong/0'; % np=6923, nice. 
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
spl(p,''); view(v);  title('guess 1'); p=setfn(p,'BCClong/b2z'); p.nc.imax=15; p.file.smod=1; 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu);
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); spl(p,''); 
huclean(p); p=resetc(p); p.sol.restart=1; p.nc.imax=5; p.sol.ds=0.01; p=pmcont(p,400);
%% BCC2tubes front from iguess and Newton 
s2=sqrt(2); zc=pi; bca=0.4; ta=1; la=0; lam=0.2; nu=1.5; % bcc2tubes, works with lz=8*lx
u=la*cos(z/s2).*(z<=zc); u=u+ta*(cos((x+y)/s2)+cos((x-y)/s2)).*(z<=zc);  
u=u+bca*((cos((x+y)/s2)+cos((x-y)/s2))+(cos((x+z)/s2)+cos((x-z)/s2))).*(z>zc); 
u=u+bca*(cos((y+z)/s2)+cos((y-z)/s2)).*(z>zc);
p.u(1:p.np)=u; p.u(p.nu+1)=lam; p.u(p.nu+2)=nu; r=norm(resi(p,p.u),'inf'); 
spl(p,''); view(v);  title('guess 2'); p=setfn(p,'BCClong/b2t'); p.nc.imax=15; 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); spl(p,''); 
%% cont of obtained solution for 1 step, then meshada, then cont further 
clf(2); p=resetc(p); p.sol.restart=1; p.nc.imax=5; p.sol.ds=0.01; p.file.smod=1; 
p=cont(p,1); op=troptions3D(); % load default trullerup-options, then overload some 
op.verbose=2; op.qualP=2.2;  op.innerit=3; op.setids=@setidsbar; 
p.sw.ips=2; p.trop=op;  p.sw.trul=1; p=oomeshada(p,'ngen',2); p=cont(p,140); 
%% a-posteriori compute stability for b2z, here just report # unstable Evals in inegli 
aux=[]; ptli=10:10:100; aux.changefile=0; inegli=fspcalc('BCClong/b2z',ptli,aux) 
%% a-posteriori compute stability for b2t, with update on disk  
ptli=1:100; aux.changefile=1; inegli=fspcalc('BCClong/b2t',ptli,aux); 
%% generate (pure) BCC and tube branches over minimal domain  
p=[]; lx=sqrt(2)*pi/2; ly=lx; lz=lx; ndim=3; lam=-0.001; nu=1.5; par=[lam; nu];  
nx=3; sw.sym=1; sw.ref=2; dir='BCClong/0s'; p.plot.pstyle=2; 
p=shinit(p,nx,lx,ly,ndim,par,lz,sw); p.pdeo.grid.plotFaces; p.np, pause 
p.sol.ds=0.001; p.sol.dsmax=0.01; p.sol.dsmin=0.001; p=setfn(p,dir); p.plot.pstyle=3; 
p.sw.bifcheck=2; p.pm.resfac=1e-3; p.usrlam=[0 0.2]; p.nc.lammax=0.4; p=cont(p,10); 
%% qswibra; 3 tubes as kernel; afterwards set some switches 
aux=[]; aux.isotol=1e-16; aux.m=3; p0=qswibra('BCClong/0s','bpt1',aux); 
p0.nc.dsmin=0.05; p0.pm.mst=4;  p0.sw.bifcheck=0; p0.pm.resfac=1e-4; p0.sw.spcalc=1; 
p0.file.smod=5; p0.nc.tol=1e-6; p0.sw.foldcheck=0; p0.sw.spcalc=1; p0.plot.pstyle=2; 
%% select BCC and cont in hot direction, use gentau for tubes;
p=seltau(p0,1,'BCClong/b1s',2); p.sol.ds=-0.02; p=pmcont(p,40); 
p=gentau(p0,1,'BCClong/t1s'); p.sol.ds=-0.02; p=pmcont(p,30); 
%% plot BD 
fnr=3; cmp=3; figure(fnr); clf; 
plotbra('BCClong/b1s','pt30',fnr,cmp,'cl','r','lab',29); 
plotbra('BCClong/t1s','pt30',fnr,cmp,'cl','b','lab',26); 
plotbra('BCClong/b2z','pt260',fnr,cmp,'cl',p2pc('r1'),'lab',[0, 200, 250]); %,'lw',4); 
plotbra('BCClong/b2t','pt100',fnr,cmp,'cl',p2pc('r2'),'lab',[0 40 80]); 
xlabel('\lambda'), ylabel('||u||'); box on; 
%% soln plots 
v=[-60 45]; 
plotsol('BCClong/b1s','pt30',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); colorbar 'off'; pause
plotsol('BCClong/t1s','pt30',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); colorbar 'off'; pause
v=[30 45]; plotsol('BCClong/b2z','pt0',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); pause 
plotsol('BCClong/b2z','pt200',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); pause 
plotsol('BCClong/b2z','pt250',1,1,3); view(v); xlabel(''); ylabel(''); zlabel(''); 
%% 
ps=3; 
for i=[0 40 80]
plotsol('BCClong/b2t',['pt' mat2str(i)],1,1,ps); view(v); xticks('auto'); yticks('auto'); 
xlabel(''); ylabel(''); zlabel(''); axis image; pause
end 
%%
fnr=3; cmp=3; figure(fnr); clf; 
plotbra('BCClong/b2t','pt100',fnr,cmp,'cl',p2pc('r2'),'lab',[0 40 80]); 
xlabel('\lambda'), ylabel('||u||'); box on; 
