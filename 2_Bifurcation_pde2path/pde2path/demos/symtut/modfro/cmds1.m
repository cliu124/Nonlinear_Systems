%% demo modfro, C1: clear and init 
close all; format compact; keep pphome; 
p=[]; lx=25; nx=50; par=[0.2 9 0]; p=modfroinit(p,lx,nx,par,'s1'); 
p.u(p.nu+1)=0.18; p1=p;   % a=par(1)=0.2->steady front, a=0.12-> osc.front via tintfreeze
%% use tintfreeze to converge to TW, and plot velocity 
nt=2e4;pmod=500;vmod=100;dt=0.05;t1=0; vel=[]; %p.u(p.nu+1)=p.u(p.nu+1)-0.02; 
[p,t1,vel]=tintfreeze(p,t1,dt,nt,pmod,vel,vmod); 
figure(7); clf; h=plot(vel(1,:),vel(2,:)); axis tight; legend('a=0.12'); 
%% switch on ext.system, cont for 1 step, then mesh-ref, then cont further 
p.nc.ilam=[1 3]; p.nc.nq=1; p.u(p.nu+3)=vel(2,end); 
p.fuha.qf=@qf; p.fuha.qfder=@qfder; p.sw.qjac=1;  
p.sol.ds=-1e-3; p=cont(p,1); p=meshada(p,'ngen',2,'sig',0.005); p=cont(p,10); 
%% swibra to Hopf, average speed <s> as aux. variable, 
aux=[]; aux.nqnew=0; % switch off steady constraint 
aux.nqh=1;   % switch on 1 'Hopf' constraint
aux.qfh=@qfh; aux.qfhder=@qfhjac; % function handles to hopf contraints 
aux.tw=1e-6; % small weight of period T in arclength often useful 
aux.dlam=0; % use trivial initial predictor 
aux.tl=100; hodir='h1'; ds=0.1; % aux.tl=80; % use 2nd version for quick results
p=hoswibra('s1','hpt1',ds,4,hodir,aux); % branch switching 
p.hopf.flcheck=1; p.hopf.fltol=1e-2; % set some additional controls 
p.file.smod=1; p.nc.dsmax=ds; p.sw.verb=0; 
p.hopf.auxp=@speedplot; % aux function for on the fly plot of s(t) 
p.nc.ilam=1; % the primary cont. parameter 
p.hopf.ilam=3; % the 2nd active parameter, here again the speed (but now average)  
p.nc.lammin=0.05; p=belon(p,1); 
%%
tic; p=cont(p,15); toc 
%% compare Hopf-orbit to tintfreeze: load some point, and freeze again 
% OK for tl=100 and a=0.147 (pt5), significant error already at at a=0.128 (pt10) 
p=loadp(hodir,'pt5','s1'); hoplot(p,1,1); speedplot(p); ylabel(''); 
vel=[]; t1=0; t2=p.hopf.T; tl=p.hopf.tl; p.u(1:p.nu)=p.hopf.y(1:p.nu,1); 
nt=1000; dt=t2/nt; pmod=round(nt/20); vmod=5; 
[p,t1,vel,uv]=tintfreezex(p,t1,dt,nt,pmod,vel,vmod); 
figure(10); hold on; plot(vel(1,:),vel(2,:),'r'); axis tight; 
set(gca,'FontSize',p.plot.fs); xlabel('t'); legend('s','s_{freeze}'); 
%freezeplotu(p,uv,1,2); % to also plot u(t,x) from freezing 
%% Floquet-plots: 
[muv1,muv2,ind,mo,h]=floqap(hodir,'pt5'); 
%% branch plotting
bf=3; figure(bf); clf; plotbra('s1','pt10',bf,3,'cl','k');
plotbra('h1',bf,3,'cl','r','lab',[5 10]); 
xlabel('a'); ylabel('s'); 
%axis([0.085 0.18 -0.29 -0.26]); 
%% soln plots
dir='h1'; pt='pt'; ind=[5 10]; clf(10); 
for i=1:length(ind);
    q=loadp(dir, [pt mat2str(ind(i))]); 
    hoplot(q, 1,1,1);  figure(1); view(0,90);  xlabel('x'); ylabel('t'); 
    title(['u_1 at ' dir '/' [pt mat2str(ind(i))]]); %shading interp; axis tight;
    hoplot2(q, 2,2,4);  figure(2); view(-20,30);  xlabel('x'); ylabel('t'); 
    title(['u_2 at ' dir '/' [pt mat2str(ind(i))]]); axis tight;
    speedplot(q); title(['s at ' dir '/' [pt mat2str(ind(i))]]); ylabel(''); pause
end 