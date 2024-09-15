%% genuine bdcurve, first with k=2, 5 initial steps 
nx=15; al=0; h0=0; v0=0; a0=0; k=2; par=[h0; v0; a0; al; k]; 
p=bdcurveinit(nx,par); p=setfn(p,'d2'); p.nc.ilam=4; p.bcsw=2; p=cont(p,5); 
%% some boundary refinement and coarsening, trial and error to choose sig 
p=loadp('d2','pt5');   % reload point (easier for trial and error) 
sigr=0.1; sigc=0.1; p=refineX(p,sigr); p=cont(p,2); p=degcoarsenX(p,sigc); 
%% continuation alternating with moveX, and refine and coarsen, parameters:
nis=15; ncs=1; % #inner steps (before ref/coars), #cont-steps (before more)
dt=0.1; nit=5; % stepsize and iterations in moveX
for i=1:3;     % outer loop, 
  for j=1:nis;  % inner loop, alternate moveX and cont  
      p=moveX(p,dt,nit); pplot(p,20); p=cont(p,ncs); 
  end 
  p=refineX(p,sigr); p=degcoarsenX(p,sigc); % refine and coarsen
end
%% k=3; 
par=[h0; v0; a0; al; 3]; p=bdcurveinit(nx,par); p.bcsw=2;  
p=setfn(p,'d3'); p.nc.ilam=4; p=cont(p,5); 
%%
p=loadp('d3','pt5'); sig=0.1; 
p=refineX(p,sig); p=cont(p,2); p=degcoarsenX(p,0.1); 
%%
for i=1:4; 
  for j=1:10; p=moveX(p,0.1,5); pplot(p,20); p=cont(p,1); end 
  p=refineX(p,sig);p=degcoarsenX(p,0.1,5); 
end
%% branch plot 
mclf(3); plotbra('d2','pt51',3,7,'lab',[30,50]); plotbra('d3','pt45','lab',[45],'cl','b'); 
xlabel('\alpha'); ylabel('A'); 
%% soln plots 
p2pglob.vi=[40 25]; p2pglob.tsw=0; p2pglob.cb=0;  p2pglob.cm='parula'; p2pglob.showbd=2; 
plotsol('d2','pt30'); pause; plotsol('d2','pt50'); pause 
plotsol('d3','pt45'); pause; p2pglob.cb=1;  plotHK('d3','pt45'); 
%% MCF test 
p=loadp('d2','pt20','dummy'); pXf=1+0.2*(rand(p.np,1)-0.5); pXf(p.idx)=1; 
p.X(:,3)=pXf.*p.X(:,3); pplot(p,10); p.fuha.flowf=@mcff; 
t=0; ts=[]; dt=0.001; ns=200; nplot=20; [p.X,t,ts]=geomflow(p,t,ts,dt,ns,nplot); 
%% show H, generally close to 0, except possibly at boundaries 
p2pglob.cb=1; plotHK(p); 
%% 
p=loadp('d2','pt50'); p.fuha.sGjac=@sGjac; jaccheck(p,0.5); 