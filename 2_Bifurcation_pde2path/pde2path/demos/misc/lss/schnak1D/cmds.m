%% 1 - initialising the problem
close all; keep pphome; p=[]; par=[sqrt(60)*sqrt(3-sqrt(8))+0.2, -0.6, 60]; 
nper=80; npp=60; p=schnakinit(p,nper,npp,par); % also use nper=80, 120, 160, 200
%% 2 - continue trivial branch to find BP 
tic; p=findbif(p,1); toc 
%% 3 - switch to periodic branch and continue. For comparison of \ and 
% lssbel, switch off stuff not related to lss 
p=swibra('p','bpt1','b1',0.02); p.sw.spcalc=0; p.sw.foldcheck=0; p.sw.bifcheck=0; 
p.sw.verb=2; p0=p; t1=tic; p=cont(p,50); t1=toc(t1); % cont with default settings  
p=p0; bw=0; beltol=1e-4; belmaxit=5; p=setbel(p,bw,beltol,belmaxit,@lss); % lssbel 
t2=tic; p=cont(p,50); t2=toc(t2); 
fprintf('t1=%g, t2=%g\n', t1,t2); plotsol(p,1,1,1); 
%% 4 - plots
figure(3); clf; plotbra('p',3,0); plotbra('b1','pt32',3,0,'lab',[21,32],'cl','b'); 
xlabel('\lambda'); plotsol('b1','pt32'); 
x=[32 64 96 128 160]; y1=[3 11 22 38 69]; y2=[1.1 2 2.8 3.9 4.8];
plot(1000*x,y1,'*-',1000*x,10*y2,'r*-'); legend('lss','10*lssbel'); xlabel('n_u'); 