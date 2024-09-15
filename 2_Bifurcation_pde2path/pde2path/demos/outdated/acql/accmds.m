% command templates for qlAC; run cell-by-cell; 
close all; keep pphome; 
p=[]; nx=40; p=acinit(p,nx); p.sw.sfem=0; screenlayout(p); p=setfn(p,'p'); 
p.sw.jac=1; p=findbif(p,1); % find bifpoints from trivial branch: 
%% one bifurcating branch 
p=swibra('p','bpt1','b1',-0.1); p.nc.dsmin=1e-5; 
p.sw.jac=1; tic; p=cont(p,10); toc 
%% works, but not checked in detail 
p=meshada(p); 

