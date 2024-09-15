% command templates for TwoFluid demo
% Publication: D. Zhelyasov, D. Han-Kwan, J.D.M. Rademacher, Global
% stability and local bifurcations in a two-fluid model for tokamak plasma,
% Preprint 2014
% run cell-by-cell; 
close all; keep pphome; 
p=[]; p=tf_init(p,25,25,2,2); 
%% find initial guess for solution via tint:
x=p.mesh.p(1,1:p.nu/3)'; y=p.mesh.p(2,1:p.nu/3)'; 
p.u(1:p.nu)=0.05*[cos(x).*sin(2*y);cos(x).*sin(2*y);cos(x).*sin(2*y)];
p=tint(p,0.1,1000,100);
%% correct it to a solution:
[p.u,res]=nloop(p,p.u); res, plotsol(p,1,3,1); 
%% for illustration: attempt continuation without phase condition:
q=p; q.sw.bifcheck=0; q.sol.ds=-0.04;  q=cont(q);
%%
clf(2); p=loadp('p','pt0'); p.sw.bifcheck=1;
%% add phase condition and continue:
p.nc.nq=1; p.fuha.qf=@tf_qf; p.fuha.qfder=@tf_qjac; 
p.sol.restart=1; p.nc.ilam=[1;5]; p.sol.ds=-1e-3; p=cont(p,20); 
%% Switch branches
p=swibra('p','bpt1','b1'); p.sw.bifcheck=0; 
p.nc.dsmax=0.05; p.sw.spcalc=1;
%%
p.sw.qjac=1;p=cont(p,10);
%%
p=swibra('p','bpt2','b2'); p.sw.bifcheck=0; 
p.nc.dsmax=0.05; p.sw.spcalc=1; p=cont(p,10);
p=swibra('p','bpt3','b3'); p.sw.bifcheck=0; p.sw.qjac=0;
p.nc.dsmax=0.05; p.sw.spcalc=1; p=cont(p,10);
%% plot BD, L^2 over delta
figure(3); clf; cmp=0;clf; ms=5;
plotbraf('p',3,cmp,'ms',ms);
plotbraf('b1','pt10',3,cmp,'cl','b','lab',5);
plotbraf('b2',3,cmp,'cl','g');
plotbraf('b3',3,cmp,'cl','r');
xlabel('\delta'); axis([0.1 0.2 0 0.5]);
%% plot BD, s over delta 
figure(3); clf; cmp=5;clf; ms=5;
plotbraf('p',3,cmp,'ms',ms);
plotbraf('b1','pt10',3,cmp,'cl','b','lab',5);
plotbraf('b2',3,cmp,'cl','g');
plotbraf('b3',3,cmp,'cl','r');
xlabel('\delta'); ylabel('s'); axis([0.1 0.2 -0.1 -0.04]);
%% plot some solutions
plotsolf('p','bpt1',1,1,2);xlabel('');ylabel('');pause
%plotsolf('p','bpt2',1,1,2);xlabel('');ylabel('');pause
plotsolf('b1','pt5',1,1,2);xlabel('');ylabel('');pause
%plotsolf('b1','pt10',1,1,2);xlabel('');ylabel('');

