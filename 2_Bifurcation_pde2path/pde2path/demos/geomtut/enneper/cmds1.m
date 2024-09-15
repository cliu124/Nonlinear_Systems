%% Enneper
global p2pglob; p2pglob.showbd=2; p2pglob.cm='parula'; p2pglob.vi=[-10,25];p2pglob.cut=0;  
nx=15; h0=0; v0=0; a0=0; alpha=0.7; k=2; par=[h0; v0; a0; alpha; k];  
p=enninit(nx,par); p.nc.lammax=2; p=setfn(p,'e1'); pplot(p,1); p.n0=p.np; 
p.nc.ilam=4; p.nc.tol=1e-6; ds=0.01; p.sol.ds=ds; p.nc.dsmax=10*ds; 
p.nc.mu2=0.1; p.nc.foldtol=0.1; p.nc.neig=20; p.nc.eigref=-0.5; % pars for spec, 
% rather many evals and negative eigref to get correct stability 
p=cont(p,1); % p.sw.jac=1; p.fuha.sGjac=@sGjac; % also works 
%% alternate moveX, cont, and refine and degcoarsen, parameters:
%p=loadp('e1','pt1','du');  ds=0.01; p.sol.ds=ds; p.nc.dsmax=5*ds; 
% hit return if no convergence; the next degcoarsenX should repair things 
nis=2; ncs=1; dt=0.05; nit=1; sigr=0.3; sigc=0.1;  cit=10; p.sw.bifcheck=1; 
for i=1:6;     % outer loop, 
  for j=1:nis;  % inner loop, alternate moveX and cont  
      p=moveX(p,dt,nit); pplot(p,20); p=cont(p,ncs); 
  end  
  tho=p.th; idold=p.idx; p=refineX(p,sigr); p=thinterpol(p,idold,tho);   
  p=degcoarsenX(p,sigc,cit,1); 
end
%% swibra to lower area
p=swibra('e1','bpt1','e2',0.05); p=cont(p,2); 
%% continuation alternating with moveX, and refine and coarsen, parameters:
nis=1; ncs=2; dt=0.025; nit=1; sigr=0.4; sigc=0.1; cit=5; % nis=3->lam=1.52
for i=1:3;     % outer loop, 
  for j=1:nis;  % inner loop, alternate moveX and cont  
      p=moveX(p,dt,nit); pplot(p,20); p=cont(p,ncs); 
  end
  tho=p.th; idold=p.idx; p=refineX(p,sigr); p=thinterpol(p,idold,tho); 
  p=degcoarsenX(p,sigc,cit,1);   
end
%% branch plot 
mclf(3); c=7; ylab='A';  
plotbra('e1','pt36',3,c,'lab',[10 16 23 36]); 
plotbra('e2','pt30',3,c,'lab',[15 30],'cl','b','ms',0); 
xlabel('\alpha'); ylabel(ylab); 
%% 
mclf(4); c=6; ylab='V';
plotbra('e1','pt16',4,c,'lab',[10]); 
plotbra('e2','pt30',4,c,'lab',[15 30],'cl','b','ms',0); 
xlabel('\alpha'); ylabel(ylab); 
%% soln plots 
p2pglob.tsw=10; p2pglob.cb=0; p2pglob.showbd=2; 
p2pglob.cm='parula';  p2pglob.vi=[-10 25];
plotsol('e1','pt10'); pause; plotsol('e1','pt16');pause; plotsol('e1','pt23'); pause; 
plotsol('e1','pt36'); pause; plotsol('e2','pt15');pause; plotsol('e2','pt30');pause; 
%%
p2pglob.cb=1; p2pglob.cm='spring'; plotHK('e1','pt36'); 
%% jaccheck 
p=loadp('e1','pt30'); pplot(p,1); p.fuha.sGjac=@sGjacenn; [Gu,Gn]=jaccheck(p,0.5);
%% cont in H 
p=swiparf('e1','pt20','e1H',1); p=cont(p,20); 