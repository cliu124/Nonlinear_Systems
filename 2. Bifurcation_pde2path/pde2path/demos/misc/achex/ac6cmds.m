%% Allen-Cahn on hexagon with x-dependent Dirichlet bc
% 
%% Initial branch
close all; %clear all; 
p=[]; p=ac6init(p); p.nc.nsteps=40; p.nc.dsmax=0.2; p=cont(p); 
%% Switch to bifurcating branches
p=swibra('p','bpt1','q',-0.1); p=cont(p);
p=swibra('p','bpt2','r1',-0.1); p=cont(p);
%% Plot bif-diagram
figure(3);clf; cmp=0; 
plotbraf('p',3,cmp,'cl','b');
plotbraf('q',3,cmp,'lab',10,'labo',0,0.2);
plotbraf('r1',3,cmp,'cl','r');
xlabel('\lambda');ylabel('||u||_2');
%% Plot some solutions
plotsolf('p','pt1',7,1,4);
plotsolf('p','pt15',8,1,1); view(-20,60);
plotsolf('q','pt10',9,1,1);view(-20,60);
plotsolf('r1','pt5',10,1,1);view(-20,60);
