%% qlAC 3D, this is expensive since running on coarse mesh requires numjac
p=[]; lx=2*pi; nx=15; dim=3; par=[0.25; -0.15; 1; -0.3; 0.2]; % c0, lam, ga, del, epsi
p=acinit(p,lx,nx,par,dim);  p.np, p.nc.tol=1e-6; p=setfn(p,'3D'); p.nc.dsmax=0.05; 
p.nc.bisecmax=5; % fewer bisections for speed;  
p.sw.jac=0; %  here jac-assembly not good enough! 
p=cont(p,10);  
%% both directions at BP1 
p=swibra('3D','bpt1','3D1c',0.1); p.sw.bifcheck=0; tic; p=cont(p,15); toc
p=swibra('3D','bpt1','3D1b',-0.1); p.sw.bifcheck=0; p=cont(p,15); 
%%
p=swibra('3D','bpt2','3D2',0.1); p.sw.bifcheck=0; p=cont(p,15); 
p=swibra('3D','bpt3','3D3',0.1); p.sw.bifcheck=0; p=cont(p,15); 
%% soln plot 
plotsol('3D1a','pt15',1,1,1); axis image; pause; 
plotsol('3D1b','pt15',1,1,1);  axis image; pause; 
plotsol('3D3','pt15',1,1,2,'alpha',0.5);  axis image; 
%% BD plot
figure(3); clf; pcmp=0; plotbra('3D','bpt3',3,pcmp,'cl','b'); 
plotbra('3D1a','pt15',3,pcmp,'cl','k','lab',15); 
plotbra('3D1b','pt15',3,pcmp,'cl','k','lab',15); 
plotbra('3D2','pt15',3,pcmp,'cl','r','lsw',0); 
plotbra('3D3','pt15',3,pcmp,'cl','m','lab',[15]); 
xlabel('\lambda'); ylabel('||u||_2'); axis([-0.1 0.21 0 19]); 
%%
p=swibra('3D','bpt1','3D1c',0.1); p.sw.jac=1; p.jacsw=1;    p.sw.bifcheck=0; tic; p=cont(p,15); toc
%% jaccheck; rel error is about 0.0002, but cont doesn't work  
p=loadp('3D1a','pt15'); p.jacsw=1; 
[Gu,Gn]=jaccheck(p); Gd=abs(Gu-Gn); e1=max(max(Gd)); spy(Gd>e1/2); 