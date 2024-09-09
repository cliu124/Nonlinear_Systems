%% sh35 on (half) disk; init, and saving of u(eps=0) as 0h/pt1 
% uses quadratic FEM, hence first run 'sethofem.m' in demos/hofem 
p=[]; dswitch=2; lam=-0.01; nu=1.4; q=1; par=[lam nu q 0 0 0]; rad=14; nref=5; 
p=shinit(p,dswitch,par,rad,nref); p=setfn(p,'0h'); % init, and set out-dir 
p.sol.ds=0.1; p=cont(p,1); % one dummy-step (to eps=0) 
%% bif directions 
aux.m=6; aux.besw=0; p0=cswibra('0h','pt1',aux);  p0.nc.eigref=-1;
%% radial 
p=gentau(p0,[0 0 0 0 1],'r'); p.sw.bifcheck=0; p.u(p.nu+2)=2; 
p.sol.ds=-0.01; p.nc.dsmax=0.02; p=cont(p,5); p=qyon(p); p=cont(p,75);  
p.nc.dsmax=0.04; p=cont(p,50);  % increase max stepsize and cont further
%% D_4^- mode 
p=gentau(p0,[0 0 0 0 0 1],'d4m'); p.sw.bifcheck=0; p.u(p.nu+2)=2; 
p.sol.ds=-0.01; p.nc.dsmax=0.02; p=cont(p,5); p=qyon(p); p=cont(p,50);  
p.nc.dsmax=0.04; p=cont(p,40); 
%% branch plot 
f=11; figure(f); clf; plotbra('r','pt130','cmp',11,'cl','k','lab',[10 120],'wnr',f,'fms',0,'ms',0); 
plotbra('d4m','pt140','cmp',11,'cl','r','lab',[20 60 90 130],'wnr',f,'fms',0,'ms',0); ylabel('||u||_2'); 
%% soln plot 
plotsol('r','pt10'); solstyle('10'); pause; plotsol('r','pt120'); solstyle('120'); pause 
plotsol('d4m','pt20'); solstyle('20','r'); pause; plotsol('d4m','pt60'); solstyle('60','r'); pause; 
plotsol('d4m','pt90'); solstyle('90','r'); pause; plotsol('d4m','pt130'); solstyle('130','r'); 
%% -------  wall mode (daisy) with nu=1.4; initially with bifdetec, then without ------------
p=gentau(p0,[0 0 1],'wall'); p.nc.eigref=-1; p.sw.bifcheck=1; 
p=cont(p,10); p.sw.bifcheck=0; p=cont(p,20);
%% daisysnake 
clf(2); p=swibra('wall','bpt1','dsnake'); p.nc.dsmax=0.01; p.sw.bifcheck=0; p=cont(p,200);
%% branch plot 
f=11; figure(f); clf; plotbra('wall','cmp',11,'cl','b','lab',[20],'wnr',f,'lp',20,'fms',0,'ms',0); 
plotbra('dsnake','cmp',11,'cl',p2pc('g1'),'fplab',[4 12 20],'wnr',f); ylabel('||u||_2'); 
%% soln plots 
plotsol('wall','pt20'); solstyle('20'); pause; plotsol('dsnake','fpt4'); solstyle('FP4'); pause 
plotsol('dsnake','fpt12'); solstyle('FP12'); pause; plotsol('dsnake','fpt20'); solstyle('FP20'); 
%% -------- FP-continuation to illustrate daisy-snake breakup at nu=1.5 ---------------
figure(2); clf; p=spcontini('dsnake','fpt12',2,'f12'); p.sol.ds=0.01;  p.sol.dsmax=0.01; 
p.plot.bpcmp=1; p.sw.bifcheck=0; p.usrlam=1.5:0.1:2; p.file.smod=1; p.sw.para=2; 
p.fuha.spjac=@spjac; p.fuha.lss=@lss; p.fuha.blss=@blss; p=cont(p,20); 
%% FP-cont plots  
f=7; c=1; figure(f); clf; plotbra('f12','pt20',f,c,'cl','m', 'lab',[12]); 
xlabel('\nu'); ylabel('\epsilon'); plotsol('f12','pt12'); nolti;
%% return to eps-cont 
p=spcontexit('f12','pt12','w1b'); p.sol.ds=1e-2; p.nc.tol=1e-4; clf(2); p.sw.spcalc=1;
p.plot.bpcmp=11;  p.file.smod=5; p.sw.foldcheck=0;  p=cont(p,100); 
%% other direction 
p=loadp('w1b','pt5','w1br'); p.sol.ds=-p.sol.ds; p=cont(p,100); 
%% branch plots  
f=11; c=11; figure(f); clf; plotbra('w1b','pt100',f,c,'cl',p2pc('g1'), 'lab',50); 
plotbra('w1br','pt100',f,c,'cl',p2pc('r1'), 'lab',100); 
xlabel('\epsilon'); ylabel('||u||_2'); 
%% soln plots 
plotsol('w1b','pt50'); solstyle('50'); pause; plotsol('w1br','pt100'); solstyle('100'); 
%% ---------- Solutions on full disk: create structure for full disk, then mirror
% half-disk soln to full disk, switch on pertinent phase conditions, then cont 
p0=[]; dsw=1; rad=14; nref=5; p0=shinit(p0,dsw,par,rad,nref);
%% daisies, need only rotational PC 
p=h2fdisk('dsnake','pt30',p0,'dsnakefull'); 
p=resetc(p); p.sw.bifcheck=0; p=qroton(p); p.sol.ds=-p.sol.ds; p=cont(p,30); 
%% radial, needs x and y phase-conditions: 
p=h2fdisk('r','pt10',p0,'rfull'); 
p=resetc(p); p.sw.bifcheck=0; p=qxyon(p); p.sol.ds=-0.02; p.nc.dsmax=0.05; 
p.nc.eigref=-1; p=cont(p,20); 
%% D_4^-,  needs all 3 phase-conditions: 
p=h2fdisk('d4m','pt20',p0,'d4mfull'); 
p=resetc(p); p.sw.bifcheck=0; p=qxyron(p); p.sol.ds=-0.02; p.nc.dsmax=0.05; 
p.nc.eigref=-1; p=cont(p,20); 
%% compare half-disk and full-disk branches
figure(3), clf; 
plotbra('dsnake','cmp',11); plotbra('dsnakefull','cmp',11,'cl','m');  ylabel('||u||_2'); 