%% cont in alpha, BCs X_3-alpha*sin k\phi=0; MUST fail when X gets 'vertical' 
nx=15; h0=0; v0=0; a0=0; al=0; k=2; par=[h0; v0; a0; al; k]; 
p=bdcurveinit(nx,par);  p=setfn(p,'b2'); p.nc.ilam=4;  p.bcsw=1; ds=0.05; 
p.sol.ds=ds; p.nc.dsmax=ds; p.nc.lammax=1.5; p.nc.lammin=-1.5;  p=cont(p,13); 
%% refine only at boundary for illustration, ... doesn't help much 
p=loadp('b2','pt6','b2r'); p.fuha.e2rs=@e2rsbdry; p.sw.rlong=1; 
p=refineX(p,sig); p=resetc(p); p=cont(p,8); 
%% like first cell, but with different ds, giving a different continuation 
nx=15; h0=0; v0=0; a0=0; al=0; k=2; par=[h0; v0; a0; al; k]; 
p=bdcurveinit(nx,par);  p=setfn(p,'b2b'); p.nc.ilam=4;  p.bcsw=1; ds=0.1; 
p.sol.ds=ds; p.nc.dsmax=ds; p.nc.lammax=1.5; p.nc.lammin=-1.5;  p=cont(p,13); 
%% back to cont in H, but BCs unchanged, e2rs back to area
p=swiparf('b2','pt6','b2H',1); p.sol.ds=0.1; p.nc.dsmax=0.1; p.file.smod=5;
p=cont(p,15); p.fuha.e2rs=@e2rsA; p=refineX(p,0.2); p=cont(p,5); 
%% other direction 
p=swiparf('b2','pt6','b2Hb',1); p.sol.ds=-0.1; p.nc.dsmax=0.1; p.file.smod=5;
p=cont(p,15); p.fuha.e2rs=@e2rsA; p=refineX(p,0.2); p=cont(p,5); 
%% branch plot of A over alpha, and some solution plots 
c=[4 7]; mclf(5); plotbra('b2','pt12',5,c,'cl','k','lab',[9 12]); 
plotbra('b2b','pt6',5,c,'cl','r','lab',5); 
plotbra('b2r','pt6',5,c,'cl','b','lab',5); xlabel('\alpha'); ylabel('A');
%% soln plot 
p2pglob.vi=[20, 40]; p2pglob.tsw=0; p2pglob.edc='k'; 
plotsol('b2','pt9'); pause; plotsol('b2b','pt5'); pause 
plotsol('b2','pt12'); axis([0 1 -1 1 -0.5 0.5]); pause; 
plotsol('b2r','pt5'); view(0,90); axis([0 1 0 1 -0.5 0.5]);
%% branch plot A over H 
c=[1 7]; mclf(5); plotbra('b2H','pt20',5,c,'cl','k','lab',[10 20]); 
plotbra('b2Hb','pt20',5,c,'cl',p2pc('gr1'),'lab',[20]); xlabel('H'); ylabel('A');
%% soln plots 
plotsol('b2H','pt10'); pause; plotsol('b2H','pt20'); pause; plotsol('b2Hb','pt20'); 
%% like first cell, but with different ds, giving a different continuation 
nx=15; h0=0; v0=0; a0=0; al=0; k=2; par=[h0; v0; a0; al; k]; 
p=bdcurveinit(nx,par);  p=setfn(p,'b2b'); p.nc.ilam=4;  p.bcsw=1; ds=0.1; 
p.sol.ds=ds; p.nc.dsmax=ds; p.nc.lammax=1.5; p.nc.lammin=-1.5;  p=cont(p,13); 
%% branch plot of A over alpha, and some solution plots 
c=[4 7]; mclf(5); plotbra('b2','pt12',5,c,'cl','k','lab',[6 12]); 
plotbra('b2b','pt10',5,c,'cl','r','lab',[5 10]); 
%% k=1 and k=4 
nx=20; h0=0; v0=0; a0=0; al=0; k=1; par=[h0; v0; a0; al; k]; 
p=bdcurveinit(nx,par); p=setfn(p,'b1'); p.bcsw=1; p=cont(p,20); pause 
k=4; par(5)=k; p=bdcurveinit(nx,par); p=setfn(p,'b4'); p.bcsw=1; p=cont(p,8);
%% branch plot of A over alpha, and some solution plots 
c=[4 7]; mclf(5); plotbra('b1','pt20',5,c,'cl','b','lab',19); 
plotbra('b4','pt7',5,c,'cl',p2pc('gr1'),'lab',7); xlabel('\alpha'); ylabel('A');
%% 
p2pglob.showbd=2; p2pglob.edc='none'; p2pglob.vi=[60,20]; p2pglob.cb=0; 
plotsol('b4','pt7'); pause; plotsol('b1','pt19'); 