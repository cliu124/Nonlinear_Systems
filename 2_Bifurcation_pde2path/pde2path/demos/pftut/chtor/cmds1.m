close all; keep pphome; 
%% CH on torus 
dir='h'; nx=30; 
dsmin=0.01; dsmax=0.1; eps=0.2; m=-0.8; lam=0; 
p=[]; R=0.5; rho=0.25; par=[m eps lam R rho 0]; % mass, eps, lagr for mass, R, rho, rot-speed 
lx=pi; ly=pi; p=chtorinit(p,lx,ly,nx,par); huclean(p); 
p.sw.verb=2; p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); %p.nc.mu2=0.01; 
p.file.smod=5; p.nc.ilam=[1 3]; p0=p; 
%% cont hom branch, with mass constraint 
p=p0; p.nc.mu1=1; p.nc.nq=1; p=findbif(p,6); p=cont(p,10); 
%% 2 vertical rings (double, but trivial) 
p=swibra(dir,'bpt1','b1',-0.05); p.nc.dsmax=0.1; p.usrlam=[-0.5 0]; p=cont(p,2); torplot(p); pause; 
p.nc.nq=2; p.nc.ilam=[1 3 6]; p=cont(p,30);  torplot(p);
%% 4 vertical rings (double, but trivial) 
p=swibra(dir,'bpt2','b2',0.05); p.nc.dsmax=0.1; p=cont(p,5); torplot(p); 
p.nc.nq=2; p.nc.ilam=[1 3 6]; p=cont(p,30);  torplot(p);
%% 2 horizontal rings (in-out)
p=swibra(dir,'bpt3','b3',0.05); p.nc.dsmax=0.05; p=cont(p,40); torplot(p);
%% 2 horizontal rings (top-bottom) 
p=swibra(dir,'bpt4','b4',0.05); p.nc.dsmax=0.05; p=cont(p,40); torplot(p);
%% spots
p=swibra(dir,'bpt5','b5',0.05); p.nc.dsmax=0.1;  p=cont(p,5); torplot(p); pause; 
p.nc.nq=2; p.nc.ilam=[1 3 6]; p=cont(p,30);  torplot(p);
%%
figure(3); clf; plotbra('h',3,7,'cl','k','lsw',0);
plotbra('b1',3,7,'cl',p2pc('b1'),'lab',[10 20]); plotbra('b2',3,7,'cl',p2pc('b3'));
plotbra('b3',3,7,'cl',p2pc('o1')); plotbra('b4',3,7,'cl',p2pc('o3'));
plotbra('b5','pt10',3,7,'cl','m'); axis([-0.75 0.01 2.2 6.7]); ylabel('E_\epsilon'); 
%%
p=loadp('b1','pt10'); plotsol(p); torplot(p); pause; p=loadp('b1','pt20'); plotsol(p); torplot(p); pause; 
p=loadp('b2','pt14'); plotsol(p); torplot(p); pause 
p=loadp('b3','pt22'); plotsol(p); torplot(p); pause; p=loadp('b4','pt21'); plotsol(p); torplot(p); pause 
p=loadp('b5','pt8'); plotsol(p); torplot(p); 
%% meshada and cont in eps, phi-PC 
dirl={'b1', 'b2', 'b5'}; ptl={'pt20' 'pt14' 'pt8'}; ll=3; 
for i=3:ll;  dir=dirl{i}; pt=ptl{i}; 
 p=swiparf(dir,pt,[dir(1:2) 'e'],[2 3 6]); p.sol.ds=-0.01; p.fuha.e2rs=@e2rs; p.nc.ngen=3; 
 p.nc.sig=0.2; p.nc.bddisty=0.1; p.nc.bddistx=0.1; p.nc.maxt=20000; 
 p=oomeshada(p); [po,t,e]=getpte(p); p.nc.lammax=1; p.nc.lammin=0.075; 
 p=cont(p,20); torplot(p); pause 
end
%% meshada and cont in eps, no PC 
dirl={'b3', 'b4'}; ptl={'pt22' 'pt21'}; ll=2; 
for i=2:ll;  dir=dirl{i}; pt=ptl{i};  
 p=swiparf(dir,pt,[dir(1:2) 'e'],[2 3]); p.sol.ds=-0.01; p.fuha.e2rs=@e2rs; p.nc.ngen=3; 
 p.nc.sig=0.2; p.nc.bddisty=0.1; p.nc.bddistx=0.1; p.nc.maxt=20000; 
 p=oomeshada(p); [po,t,e]=getpte(p); p.nc.lammax=1; p.nc.lammin=0.075; 
 p.nc.foldtol=0.1 ; p=cont(p,20); torplot(p); pause 
end
%%
figure(3); clf; plotbra('b1e','pt8',3,7,'cl',p2pc('b1'),'lab',7); grid on; 
xlabel('\epsilon');  ylabel('E_\epsilon'); 
%%
figure(3); clf; 
plotbra('b2e',3,7,'cl',p2pc('b3'),'ms',0); %,'lab',10);
plotbra('b3e',3,7,'cl',p2pc('o1'),'lab',9); plotbra('b4e',3,7,'cl',p2pc('o3'),'lab',12);
plotbra('b5e','pt12',3,7,'cl','m','lab',11); grid on;  xlabel('\epsilon');  ylabel('E_\epsilon'); 
%axis([0.07 0.2 3.1 9]); 
%%
p=loadp('b1e','pt7'); plotsol(p); torplot(p); pause 
p=loadp('b3e','pt9'); plotsol(p); torplot(p);  pause 
p=loadp('b4e','pt12'); plotsol(p); torplot(p);  pause 
p=loadp('b5e','pt11'); plotsol(p); torplot(p); 

