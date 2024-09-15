%% MCF, with initial large H; to handle meshing, alternate flow and coarsening 
p2pglob.cut=4; p2pglob.vi=[30,40]; p2pglob.cm='spring'; % graphics settings 
p=loadp('b2Hb','pt30','dummy'); p.t=0; pplot(p,10); % load and plot initial cond.
p.fuha.flowf=@mcff; t=0; ts=[]; dt=0.001; ns=500; nplot=50; % prepare MCF 
[p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot);  pause; % flow 
p.sw.rlong=1; sig=0.3; p.fuha.e2rs=@e2rsAi; p=coarsenX(p,sig); pause % coarsen
[p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); p=coarsenX(p,sig); pause % again 
%%
[p.X,t,ts]=geomflow(p,t,ts,dt,2*ns+20,nplot); % flow further 
%% plot time series, go back here after different flows tested below
mclf(11); plot(ts(1,:),ts(2,:),ts(1,:), ts(3,:)); legend('A','V'); 
set(gca,'fontsize',12); xlabel('t'); 
%% nice and stable! 
p=loadp('b2','pt10','dummy');  pplot(p,10); p.fuha.flowf=@mcff; 
t=0; ts=[];dt=0.001; ns=200; nplot=20; [p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
%% perturb: still OK! 
p=loadp('b2','pt10','dummy'); pXf=1+0.5*(rand(p.np,1)-0.5); pXf(p.idx)=1; 
p.X(:,3)=pXf.*p.X(:,3); pplot(p,10); p.fuha.flowf=@mcff; 
t=0; ts=[]; dt=0.001; ns=200; nplot=20; [p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
%% perturb: still OK! 
p=loadp('d2','pt20','dummy'); pXf=1+0.2*(rand(p.np,1)-0.5); pXf(p.idx)=1; 
p.X(:,3)=pXf.*p.X(:,3); pplot(p,10); p.fuha.flowf=@mcff; 
t=0; ts=[]; dt=0.001; ns=200; nplot=20; [p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
%% show H, generally close to 0, except possibly at boundaries 
p2pglob.cb=1; plotHK(p); 