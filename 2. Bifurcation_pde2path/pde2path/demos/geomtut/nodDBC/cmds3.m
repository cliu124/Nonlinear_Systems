%% One inner buckle nodoid, cont in A 
close all; keep pphome; clear classes; % modification in stanpdeo2D 
%% for convenience, control behavior of pplot by p2pglob
global p2pglob; p2pglob.vi=[10 30]; p2pglob.edc='k'; p2pglob.faceal=1; p2pglob.cut=0; 
p2pglob.tsw=2; % use (A,H) in title
%%
V0=0; A0=0; h0=1.8; a=0.7; mu=0; laml=0; % mu for rot.PC, laml=length (for getV), 
par=[h0; V0; A0; a; 0; mu; laml]; % a used for init. par(5) not used
%     1                 6 
lx=pi/2-0.02; ly=pi; ny=50; nx=70;  sym=0; p=[];
p=nodinitl(p,lx,ly,nx,ny,par,sym); p=setfn(p,'lNA'); huclean(p); pplot(p); p0=p;
%% go 
p=p0; p.nc.ilam=[3 1]; p.nc.nq=1; p.nc.tol=1e-2; p.nc.lammax=20; 
p.fuha.qf=@qfA; p.fuha.qfder=@qjacA; p.sw.qjac=1;  p.sw.jac=1; p.nc.neig=20; 
p.plot.bpcmp=1; p.sol.ds=0.1; p.nc.dsmax=0.1; p.sw.para=1; p.nc.foldtol=0.05; 
p.nc.mu1=2; p.nc.mu2=0.1; p.sf=1e3;
p.nc.bisecmax=5; p=cont(p,1); p.nc.tol=1e-8; p=cont(p,20); 
%% 1st BP, simple, breaks Z2 symmetry 
p=swibra('lNA','bpt1','lNA1',0.05); p=cont(p,10); 
%% second bifurcation point, double, use gentau to choose bif direction 
% do 2 steps without PC, hence also without bifcheck, since \pa_\phi X in kernel
aux.besw=0; aux.mu2=0.5; aux.m=2; aux.ral=1; p1=qswibra('lNA','bpt2',aux); 
p=gentau(p1,[0.5,0.5],'lNA2',2); p.sol.ds=0.1; p.sw.bifcheck=0; p=cont(p,2);
%% switch on rot constraint, rot.speed is par(6) 
p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot;
p.sol.ds=0.4; p.sw.bifcheck=1; p.file.count=p.file.count+1;  p=cont(p,10); 
%% third bifurcation point 
aux.besw=0; aux.mu2=0.5; aux.m=2; aux.ral=1; p1=qswibra('lNA','bpt3',aux); 
p=gentau(p1,[0.5,0.5],'lNA3',2); p.sol.ds=0.4; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,2);
%% switch on rot constraint 
p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot;
p.nc.tol=1e-6; p.sol.ds=0.1; p.file.count=p.file.count+1; p=cont(p,10);
%% 4'th bifurcation point 
aux.besw=0; aux.mu2=0.2; aux.m=2; aux.ral=1; p1=qswibra('lNA','bpt4',aux);
p=gentau(p1,[0.5,0.5],'lNA4',2); p.sol.ds=0.2; p.nc.tol=1e-4; p.sw.bifcheck=0; p=cont(p,2);
%% switch on rot constraint 
p.nc.nq=2; p.nc.ilam=[3 1 6]; p.fuha.qf=@qfArot; p.fuha.qfder=@qjacArot; 
p.file.count=p.file.count+1; p.nc.tol=1e-6; p.sol.ds=0.4; p=cont(p,8);
 %% branch plot
f=3; c=1; mclf(f); 
plotbra('lNA',f,c,'cl','k','lab',[2 20]);
plotbra('lNA1',f,c,'cl','b','lab',[10]); 
plotbra('lNA2',f,c,'cl','r','lab',[10]); 
plotbra('lNA3',f,c,'cl','m','lab',[10]);
plotbra('lNA4',f,c,'cl',p2pc('g1'),'lab',[10]); 
%% soln plot
p2pglob.vi=[10 40]; p2pglob.cut=0; p2pglob.showbd=0; 
plotsol('lNA','pt2'); pause; plotsol('lNA','pt20'); pause;
plotsol('lNA1','pt10'); pause; plotsol('lNA2','pt10'); pause; 
plotsol('lNA3','pt10'); pause; plotsol('lNA4','pt10'); 
%% example of degcoarsenX - cont - loop 
p=loadp('lNA2','pt3','lNA2c'); p.nc.tol=1e-6; p.sw.bifcheck=0; p2pglob.cut=1; 
sig=0.2; nit=6; keepbd=1; 
for i=1:6; p=degcoarsenX(p,sig,nit,keepbd); p=cont(p,4); end
%% coarsen via coarsufu
p=loadp('lNA2','pt3','lNA2cc'); p.nc.tol=1e-6; p.sw.bifcheck=0; 
p.fuha.ufu=@coarsufu; p.nc.sig=0.3; p.nc.errbound=100; p=cont(p,10); 
%% branch plot, mesh-distortion at pos npar+4=11; 
f=3; mclf(f); c=11; % was 13 
plotbra('lNA2','pt10',f,c,'cl','r','lab',10,'fp',1); 
plotbra('lNA2c','pt22',f,c,'cl',p2pc('r2'),'lab',[14,22],'fp',1); 
plotbra('lNA2cc','pt30',f,c,'cl','m','lab',[],'fp',4); 
ylabel('\delta_{mesh}=max(h/r)'); 
%%
p2pglob.vi=[10 20]; p2pglob.cut=1; p2pglob.showbd=2; 
plotsol('lNA2','pt10'); pause; plotsol('lNA2c','pt14');pause; plotsol('lNA2c','pt22');
%% further degcoarsenX checks: 
p=loadp('lNA2','pt16','lNA2c'); o1=meshqdat(p);  
keepbd=1; nit=6; sig=1; p2pglob.cut=1;  pplot(p,10); 
p=degcoarsenX(p,sig,nit,keepbd); o2=meshqdat(p); 
pause; p.nc.tol=1e-6; p.sw.bifcheck=0; p=cont(p,5); 
%% coarsen 3042 to 2143, nice. keepbd=1; nit>5! 
p=loadp('lNA3','pt10','lNA3c'); keepbd=1; nit=6; sig=1; p2pglob.cut=1;  pplot(p,10); 
p=degcoarsenX(p,sig,nit,keepbd); pause; p=cont(p,5); 
%% refine after coarsen, no problem, but needs p.DIR cleaned 
p=loadp('lNA3c','pt15','lNA3cr'); plotsol(p); sig=0.1; p.sw.rlong=1; p.DIR=[]; 
p=refineX(p,sig); pause; p=cont(p,4); 
%% only refine, good for 5 steps 
p=loadp('lNA3','pt10','lNA3r'); plotsol(p); sig=0.1; p.sw.rlong=1; 
p=refineX(p,sig); pause; p=cont(p,10); 
%% coarsen after refine, nice. keepbd=1 also works, but very acute Ts at bdry!  
p=loadp('lNA3r','pt15','lNA3rc'); p.branch=[]; keepbd=1; nit=4; sig=0.5; 
p=degcoarsenX(p,sig,nit,keepbd); p=cont(p,5); 
%%
p2pglob.vi=[10 30]; p2pglob.cut=1; 
p=loadp('lNA3c','pt12','dummy'); p.nt, pplot(p,10);% plotHK(p); 
p=moveX(p,0.01,3); p.nt, pplot(p,1); pause; p=cont(p,2); 



