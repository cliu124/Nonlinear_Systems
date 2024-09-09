% command templates for Gross-Pitaevski and vector GP; run cell-by-cell; 
close all; keep pphome; 
%% --------------------------------- scalar GP 
lx=5; nx=20; om=0.2; mu=0.5; ps=0; rs=0; m=1; % lag-multipl., phase and rot 
par=[om mu ps rs];  sw.sym=2; % parameters, and type of mesh 
p=[]; p=gpinit(p,lx,nx,par,sw,m); p.nc.neig=20; p=setfn(p,'p');
%% continue 
p.nc.ilam=[1 3]; p.nc.nq=1; p=cont(p,40); % fix phase and cont  
%% branch plotting, compos: 1-4=pars, 5=rmax, 6=imax, 7=imax/rmax, 8=N
figure(3); clf; plotbra('p','pt40',3,5,'lab',[1,10]);
%% solution plotting 
%p=loadp('p','pt1'); plotall1(p); 
p=loadp('p','pt10'); plotall1(p); 
%% --------------------------------- vector GP 
close all; keep pphome; 
p=[]; lx=4; nx=40; sw.sym=2; om=0.2; mu1=0.5; mu2=2; ps=0; rs=0; 
par=[om mu1 ps rs mu2]; m=1; % parameters, mu2=chem.potential of 2nd compo 
p=vgpinit(p,lx,nx,par,sw,m); p.nc.neig=20; p=setfn(p,'v'); p.nc.dsmax=0.05; 
p=oomeshada(p); plotsol(p,1,3,3); 
%% here cont without phase fix
p=cont(p,20); 
%% branch plotting, 1-5=pars, 6=r1, 7=i1, 8=i1/r1, 9-11: r2..i2/r2, 12=N1, 13=N2, 14=N1+N2
figure(3);clf; cmp=11; 
plotbra('v','pt20',3,cmp,'lab',[10,20]); xlabel('\omega'); 
%% solution plotting
p=loadp('v','pt20'); plotall2(p);  