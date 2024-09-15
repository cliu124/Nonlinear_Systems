%% AC on cylinder with top lid, i.e., 2 problems on 2 domains with coupling via bdry 
close all; keep pphome; 
%% initialization 
p=[]; R=1; s=0; par=[0.1 -0.2 1 R s]; % parameters, s=speed for PC 
lx=pi; ly=pi; p=acinit(p,lx,ly,70,par); p.ly=ly; p.sw.bifcheck=2; p.nc.mu1=0.5; p.nc.mu2=0.1; 
p.sw.verb=1; p.nc.bisecmax=15; p.sf=1e3; p.nc.neig=30; p.sw.jac=1; p.nc.del=1e-6; 
p.nc.ilam=2; p.nc.lammax=2.2; p.sol.ds=0.1; p.nc.dsmax=0.51; p=setfn(p,'tr'); 
%% cont. of trivial branch 
p=cont(p,10);
%% bif to spat.homog. branch 
p=swibra('tr','bpt1','b1',0.1); p.nc.dsmax=0.11; p=cont(p,40); 
%% secondary bif 
p=swibra('b1','bpt1','b1-1',0.1); p=cont(p,20); 
%% '2nd' BP, use besw=0 trick, i.e.: just compute kernel, then follow some selected 'pure' 
% bifurcation directions; here these are double (phi-dependent), hence use PC
aux=[]; aux.m=6; aux.besw=0; p0=cswibra('tr','bpt2',aux); p0.sw.bifcheck=0; 
p0.sw.verb=2; p0.sol.ds=0.025; p0.nc.dsmax=0.2;  p0.nc.intol=0.01; p0.nc.tol=1e-8; pause 
p=gentau(p0,[1],'b2'); p=cont(p,5); p.nc.nq=1; p.nc.ilam=[2 5]; p=cont(p,30); 
p=gentau(p0,[0 0 0 1],'b3'); p=cont(p,2); p.nc.nq=1; p.nc.ilam=[2 5]; p=cont(p,40); 
%% '3rd' (not really!) BP, proceed as above
p0=cswibra('tr','bpt3',aux); p0.sw.bifcheck=0; p0.sw.verb=2; p0.sol.ds=0.01; p0.nc.tol=1e-6; pause
p=gentau(p0,[1],'b4'); p=cont(p,2); p.nc.nq=1; p.nc.ilam=[2 5]; p=pmcont(p,80); 
p=gentau(p0,[0 0 1],'b5'); p=cont(p,2); p.nc.nq=1; p.nc.dsmax=0.1; p.nc.ilam=[2 5]; p=cont(p,50); 
%% BDs 
f=3; c=6; figure(f); clf; plotbra('tr',f,c,'cl','k','lsw',0); 
plotbra('b1',f,c,'cl','r','lsw',0); plotbra('b1-1','pt20',f,c,'cl',p2pc('r2'),'lab',5);
plotbra('b2',f,c,'cl',p2pc('b1'),'lab',15); plotbra('b3',f,c,'cl',p2pc('b3'),'lab',20); 
plotbra('b4',f,c,'cl',p2pc('o1'),'lab',30); plotbra('b5',f,c,'cl',p2pc('o3'),'lab',30); 
axis([-0.5 2.1 0 8.5]); xlabel('\lambda'); ylabel('||u_1||_2'); 
%% soln plots 
plotsol('b1-1','pt5',1,10); xticks([-0.5 0.5]); yticks([-0.5 0.5]); pause; 
plotsol('b2','pt15',1,10); xticks([-0.5 0.5]); yticks([-0.5 0.5]); colorbar off; pause; 
plotsol('b3','pt20',1,10); xticks([-0.5 0.5]); yticks([-0.5 0.5]); colorbar off; pause; 
plotsol('b4','pt30',1,10); xticks([-0.5 0.5]); yticks([-0.5 0.5]); colorbar off; pause; 
%%
plotsol('b5','pt30',1,10); xticks([-0.5 0.5]); yticks([-0.5 0.5]); colorbar off; 
%% continue in some other param, here R 
p=swiparf('b3','pt20','b3-R',[4 5]); p.sol.ds=0.05; p.nc.dsmax=0.2; clf(2); p=cont(p,20); 
%% BD of R-cont 
figure(f); clf; plotbra('b3-R',f,c,'cl',p2pc('o3')); ylabel('||u_1||_2'); 
%%
plotsol('b3-R','pt8',1,10);  colorbar off;
%% Jaccheck ...
[Ga,Gn]=jaccheck(p); Gd=(abs(Ga-Gn)); e1=max(max(Gd)); spy(Gd>e1/10); 

