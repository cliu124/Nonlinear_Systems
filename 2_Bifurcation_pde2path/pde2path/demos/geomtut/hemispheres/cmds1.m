%% VPMCF tests for CMC half sphere with NBCs in z=0 plane 
% reasonable conservation of V, and convergence to steady state, which then 
% can be continued. PCs for x,y translations as cog (center of gravity). 
% Also works without PCs, then transl. evals about 1e-2.  
close all; keep pphome; 
global p2pglob; p2pglob.tsw=1; p2pglob.vi=[20, 40]; p2pglob.edc='k'; 
p2pglob.cut=0; p2pglob.cm='parula'; 
%% init; 
nx=15; h0=0; v0=0; a0=pi; sx=0; sy=0; par=[h0; v0; a0; sx; sy]; %nx=10; % initial pars 
p=hspinit(nx,par); p=setfn(p,'hs1'); p.sol.ds=0.1; p.nc.dsmax=5; 
p.sw.jac=0; p.sw.qjac=1; p.nc.ilam=[2 1]; p.nc.nq=1; p.fuha.qf=@qfV; p.fuha.qfder=@qjacV;
p.plot.bpcmp=1; p.sw.spcalc=1; p.sw.verb=2; pplot(p); 
%% one cont.step without PC, then with PC 
p=cont(p,1); p.fuha.qf=@qf2; p.fuha.qfder=@qjac2; p.nc.nq=3; p.nc.ilam=[2 1 4 5]; p=cont(p,5); 
%% refine, then cont further 
p=loadp('hs1','pt5','hs1r'); p.sw.nobdref=0; p.sw.rlong=1; p.file.smod=5; 
sig=0.25; p=refineX(p,sig); pplot(p); pause; p=cont(p,5); pause; p=cont(p,45); 
%% branch plot of H over V and of pos.multiplier sx 
c=[2 1]; mclf(5); plotbra('hs1r','pt50',5,c,'lab',50); xlabel('V'); ylabel('H'); box on; 
mclf(7); plotbra('hs1r','pt50',7,[2,5],'cl','r'); xlabel('V'); ylabel('sx'); % multiplier
%% compare A and V as functions of r; for sphere: V=r*A/3
p=loadp('hs1r','pt30'); br=p.branch; r=br(11+4,:); V=br(11+1,:); A=br(11+2,:); 
mclf(10); plot(r,r.*A/3,'linewidth',2); hold on; plot(r,V,'*'); hold off; 
xlabel('r'); legend('r*A/3','V'); axis tight; set(gca,'fontsize',14); 
% plotbra('hs1r','pt50',3,9,'cl','r'); plotbra('hs1r','pt50',3,10,'cl','b');  % xmin,xmax 
%% some solution plots 
p2pglob.cb=0; p2pglob.cm='parula'; plotsol('hs1r','pt50'); pause; plotsol('hs1r','pt50'); 
%% prepare VPMCF: pert.some X, set parameters for VPMCF 
p2pglob.cm='spring'; p2pglob.cb=1; p2pglob.vi=[20,40]; p2pglob.showbd=2; 
p=loadp('hs1r','pt50','dummy');  clf(11); 
N=getN(p,p.X); X=p.X; th=angle(X(:,1)+1i*X(:,2)); z=X(:,3); z2=max(z); 
pXf=cos(th).*(z2-z)+0.1*(rand(p.np,1)-0.5); p.X=p.X+0.4*pXf.*N; plotHK(p); 
figure(1); title('H for initial X for VPMCF'); [xpos0,ypos0]=getpos(p) 
p.fuha.flowf=@vpmcff; ts=[]; t=0; dt=0.0005; ns=round(1/dt); nplot=round(0.01/dt);
%% go. (and repeat cell for longer DNS) 
[p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
plotHK(p); figure(1); title(['H for X(' mat2str(t,2) ')']); [xpos1,ypos1]=getpos(p) 
%% plot A and V time-series 
mclf(11); plot(ts(1,:),ts(3,:),'linewidth',2); hold on; plot(ts(1,:),ts(2,:),'linewidth',2); 
legend('V','A');xlabel('t');  set(gca,'fontsize',14); ts(3,end)/ts(3,1)
%% reset V and H, then cont again 
hm=plotHK(p); V0=getV(p,p.u); p.u(p.np+2)=V0; p.u(p.np+1)=hm; p=cont(p,1); p=cont(p,5); 
%% check cont without PC; main difference are 2 'almost 0' evals 
p=loadp('hs1','pt1','hs1b'); [xpos0,ypos0]=getpos(p), pause, p.nc.dsmax=5; p.nc.dlammax=5; 
p=cont(p,20); [xpos1,ypos1]=getpos(p) 

