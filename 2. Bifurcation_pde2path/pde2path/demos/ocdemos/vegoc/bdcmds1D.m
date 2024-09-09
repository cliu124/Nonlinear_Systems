close all; keep pphome; 
%% stat. BD for vegOC, init, FSS branch 
close all; p=[]; lx=5; nx=60; ly=sqrt(3)*lx/2; sw=1; rho=0.03; 
ndim=1; p=veginit(p,lx,ly,nx,sw,rho,ndim); 
p=setfn(p,'FSS'); screenlayout(p); 
p.nc.ilam=8; p=setlam(p,23); p.sol.ds=-2; p.nc.nsteps=50;  p0=p; 
p=findbif(p,4); % find bif.points from FSS branch 
%%  delete pt-output in dir FSS generated by findbif, to keep dir small 
% and use cont to find the desired usrlam output points
delete('FSS/pt*.mat'); p=p0; p.sol.ds=-1; p.sw.bifcheck=0; p=cont(p,50); 
%% bif from FSS (set bpt* and p* and repeat as necessary) 
p=swibra('FSS','bpt1','p1',-2); p.sw.bifcheck=0; p.sw.spcalc=0;p=cont(p,75); 
%% the 2nd and 3rd branch 
p=swibra('FSS','bpt2','p2',-2); p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p); 
p=swibra('FSS','bpt3','p3',-2); p.sw.bifcheck=0; p.sw.spcalc=0; p=cont(p,35); 
%% plot the BD, max v over R  
figure(3); clf; pcmp=1; 
plotbra('FSS',3,pcmp,'cl','k','lab',[14,18,29,45]); 
plotbra('p1','pt74',3,pcmp,'cl','b', 'lab',[11,16,34,38, 49,65]); 
plotbra('p2','pt54',3,pcmp,'cl','r', 'lab',[23 37]); 
plotbra('p3','pt37',3,pcmp,'cl','m','lab',[8]); 
text(31,400,'FSS','fontsize',16);text(30,520,'p1','fontsize',16);
text(13.5,200,'p2','fontsize',16);text(24.5,570,'p3','fontsize',16);
xlabel('R'); ylabel('max v');
%% BD, J_{c,a} over R 
figure(3); clf; pcmp=4; 
plotbra('FSS',3,pcmp,'cl','k','lab',[14,18,29]); 
plotbra('p1','pt74',3,pcmp,'cl','b', 'lab',[11,16,34,38, 49]); 
xlabel('R'); ylabel('J_{c,a}'); axis([19 30 12 29]);
%% soln plots (pauses useful for exporting plots) 
plot1Df('p1','pt11',1,10,10,10,3);pause; plot1Df('p1','pt16',1,10,10,10,3);pause 
plot1Df('p1','pt34',1,10,10,10,3);pause; plot1Df('p1','pt38',1,10,10,10,3);pause 
plot1Df('p1','pt49',1,10,10,10,3);
%% values for table 
valf('FSS','pt22'); valf('FSS','pt18'); valf('FSS','pt29'); 
valf('FSS','pt45'); 



