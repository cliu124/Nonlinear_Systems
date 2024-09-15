keep pphome; 
%The first thing to do is to create the double domain
p=[];
lx=14;%%This is important. What is the radius of the solution to mirror
ly=lx;
nx=0;
ndim=9; % ndim=9 means secpdeo2 
lam=-0.01;
nu=2; 
q=1;
par=[lam nu q 0 0 0];
sw.sym=2;
aux.rmin=0; 
aux.nr=80;
aux.al=1; % # and scaling of discr.points in r
aux.nphi=round(aux.nr); % # discr.points in phi (not used anymore) 
aux.phi=pi/2;
aux.nr=80; 
aux.jf=4; 
dir='pis'; % 
p=shinit(p,nx,lx,ly,ndim,par,sw,aux);
p=setfn(p,dir);
huclean(p);
figure(2);
grid on;
plotsol(p,1,1,0);
p.np 
p.sw.bifcheck=2;
p0=p; 
%% dummy refine near 0 
p.tau=0*p.u; p.nc.ngen=1; p.nc.maxt=40000; p.nref=600; 
%p.pdeo.grid.rlong=1;
p=oomeshada(p); 
%% do one step for saving 
p=cont(p,1);
%% now double 
p=loadp('pis','pt1');
p.al0=1; 
q=loadp('nidat','pt10'); plotsol(q,1); plotsol(p,2); 
%This is a solution with phi=pi/4
 %p1=q2p(p,q,pi/4);plotsol(p1,1) %it's the reflexion, but not what we want 
 %p=q2p(p,q,-pi/4);plotsol(p,1); %Still not what we want
 p=q2p(p,q,-pi/4);plotsol(p,11);plotsol(p,12,2); %I think this is ok
 radperi(q,31);radperi(p,32); 
 p=setfn(p,'dd1b');
 %plotsol(p); 
 p.sw.bifcheck=0; p.sw.foldcheck=0; p0=p; clf(2); % save for tests
 %% 
 p.sol.ds=-0.01; p.nc.dsmax=0.05;  p=cont(p,150);
 %%
 labli=[10 40 60 100 140]; 
 figure(3); clf;  plotbra('dd1b','fp',10,'lab',labli); 
 %% soln plots 
br='dd1b'; nl=length(labli);
for i=1:nl; p=loadp(br,['pt' mat2str(labli(i))]); 
    plotsol(p,1,1,2); axis image; nola; rplot(p,20,2*aux.nr,10); pause; 
end 
 %% 
 p=p0;  p=qron(p); 
 p=cont(p,5); 
 %% a yet smaller dsmax, same branch 
 p=loadp('dd1b','pt30','dd1c'); p.nc.dsmax=0.02; 
 p=cont(p,150); 
 %% checking pmcont, same branch ! 
 p=loadp('dd1b','pt30','dd1d'); p.nc.dsmax=0.02; 
 p=pmcont(p,350); 
