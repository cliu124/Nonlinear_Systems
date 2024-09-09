%% continuation of a curved channel in eta1 yields pearling 
p=swiparf('m1','pt10','m1e',[1 2 6]); p.sol.ds=0.01; p.nc.lammax=3; p.nc.dsmax=0.1; 
p=cont(p,50); 
%% put BPs here for inspection 
p=swibra('m1e','bpt2','du',0.01);
%% swibra and cont 
p=swibra('m1e','bpt1','m1ep1',0.01); p.nc.dsmax=0.05; p=cont(p,20);
%% plot BD, first at low m, gamma=3 
figure(3);clf;cmp=8; plotbra('m1e','pt20',3,cmp,'cl','k','bplab',1); 
plotbra('m1ep1','pt20',3,cmp,'cl','b','lab',20); 
xlabel('\eta_1'); ylabel('max(u)'); 
%% solution plots 
ps=2; myax='image'; v=[0,90]; 
plotsol('m1e','bpt1',10,1,ps); axis(myax); view(v); nolti; pause 
plotsol('m1ep1','pt20',10,1,ps); axis(myax); view(v); nolti; pause
%% tangent plot via calling swibra again, then some settings
q=swibra('m1e','bpt1','du',0.01); nolti; view(0,90); axis image; 
title('\tau_1 at ebp1','FontSize',14); set(gca,'FontSize',14);
%% time-integration from a perturbed channel (with fixed mass) 
p=loadp('s1','bpt4','t1'); m=p.u(p.nu+4), m=-0.8; % target mass
p.u(1:p.np)=p.u(1:p.np)+0.1*(rand(p.np,1)-0.5); % random pert. 
mp=sum(p.mat.Ms*p.u(1:p.np))/p.Om; p.u(1:p.np)=p.u(1:p.np)+(m-mp); % set mass 
up; plotsol(p,1,1,1); % plot sol. and prepare DNS 
t1=0; ts=[]; nc=0; dt=0.01; nt=5000; pmod=250; smod=2500; p.mat.Kadv=0; 
%% the DNS loop, repeat this cell until near steady state 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalft); axis image; nolti
%% run Newton on solution to go to 'steady state' 
p.u(p.nu+4)=m; % set 'contraint mass' to target mass
[u1,r1,i1,Gu,Glam,p]=nloop(p,p.u); fprintf('res=%g, iter=%g\n',r1,i1); 
p.u=u1; plotsol(p,1,1,p.plot.pstyle); 
p=setfn(p,'sp1a'); p=resetc(p); stansavefu(p); p.sw.verb=1; 
%% continue (if desired)
p=loadp('sp1a','bpt1'); p.sw.bifcheck=0; p.sw.spcalc=1; plotsol(p); 
p=resetc(p); p.sol.restart=1; p.sol.ds=-0.01; p=cont(p,80);