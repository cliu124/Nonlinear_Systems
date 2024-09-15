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