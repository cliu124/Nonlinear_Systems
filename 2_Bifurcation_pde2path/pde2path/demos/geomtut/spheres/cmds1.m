%% cont sphere in H and vol (trivial), then check VPMCF
% unclear: this needs a few initial steps until H=H(r) becomes correct 
% way out: start at somewhat lower r0, throw away initial points
global p2pglob; p2pglob.showN=0; p2pglob.cb=1; r0=0.25; 
h0=0; v0=0; a0=0; sx=0; sy=0; sz=0; % lagr.multipl.for pos. constraints (not used)
par=[h0; v0; a0; sx; sy; sz]; sw=3; fn=['S' mat2str(sw)]; % discretization 
p=sphereinit(par,r0,sw); p=setfn(p,fn); plotsol(p); [x,y,z]=getpos(p), pause 
p.sw.verb=2; p.sw.spcalc=1; p.nc.neig=3; p.sol.ds=0.01; p=cont(p,30); plotHK(p); 
%% branch plot of H over V, and some solution plots. 
c=[2 1]; mclf(5); plotbra(fn,'pt30',5,c,'lab',30,'fp',11); 
xlabel('V'); ylabel('H'); box on; plotHK(fn,'pt30'); 
%% compare A and V as functions of r; for sphere: V=r*A/3, fits perfectly! 
p=loadp(fn,'pt20'); br=p.branch; r=br(6+10,:); V=br(6+7,:); A=br(6+8,:); 
mclf(10); plot(r,r.*A/3,r,V,'*'); xlabel('r'); legend('r*A/3','V'); 
axis tight; set(gca,'fontsize',14); grid on; box on; 
%% VPMCF, load pt, perturb, then go 
p2pglob.cm='spring'; p2pglob.Hf=1; p2pglob.cut=0; p2pglob.showN=0; 
p=loadp(fn,'pt30'); X=p.X; th=angle(X(:,2)+1i*X(:,3)); z=X(:,3); z2=max(z); 
x=X(:,1); x2=max(x); pXf=sin(th).*(abs(x)-x2)+0.2*(rand(p.np,1)-0.5); N=getN(p,p.X); 
p.X=p.X+0.4*pXf.*N; plotHK(p); figure(1); title('H for initial X for VPMCF'); pause
p.fuha.flowf=@vpmcff; ts=[]; t=0;dt=0.005;ns=round(1/dt);nplot=round(0.05/dt);
[p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
plotHK(p); figure(1); title(['H for X(' mat2str(t,2) ')']); 
%% repeat flow (if desired) 
[p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
plotHK(p); figure(1); title(['H for X(' mat2str(t,2) ')']); 
%% plot A and V time-series 
mclf(11); plot(ts(1,:),ts(3,:),'linewidth',2); hold on; 
plot(ts(1,:),ts(2,:),'linewidth',2); legend('V','A');xlabel('t');  set(gca,'fontsize',14); 
ts(3,end)/ts(3,1)  % assess change in V 