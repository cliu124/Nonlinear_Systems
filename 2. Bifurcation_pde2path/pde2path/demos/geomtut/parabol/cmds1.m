%% paraboloid, linear, for testing different meshes and mesh-adaptations 
keep pphome; hucl; 
global p2pglob; p2pglob.edc='k'; p2pglob.cut=0; p2pglob.showN=0; p2pglob.cb=1; 
p2pglob.axlab=1; p2pglob.showbd=0; 
%% init, nx=initial mesh resol., ssw=sym-switch, nref=#ref steps for conv tests
lx=1; ly=1; a=2; b=1; ssw=[]; ssw.sym=0; nx=6; nref=4; par=[a;b;0;0]; p=[]; 
p=parabolinit(p,lx,ly,par,ssw,nx); pplot(p); p0=p; % for reusing p 
fig(1); colorbar off; nola; title([]); % initial plot; 
% uncomment next line for plot if ssw.sym=1 and subsequent retrig 
%pause; p1=retrigX(p); pplot(p1); fig(1); colorbar off; nola; title([]); 
%% strategy 0, M_full, rlong=1, just refine (e2rsA) 
p=p0; p.sw.msw=1; p.sw.rlong=1; p=solfixpar(p); ev0=eplot(p); 
for i=1:nref;  p=refineX(p); p.X=setbcX(p); 
 [p,res]=solfixpar(p); if res<p.nc.tol; ev0=[ev0 eplot(p)]; end 
end 
%% strategy 1, M_Vor, rlong=1, just refine 
p=p0; p.sw.msw=0; p.sw.rlong=1; p=solfixpar(p); ev1=eplot(p); 
%
for i=1:nref;  p=refineX(p); p.X=setbcX(p); 
 [p,res]=solfixpar(p); if res<p.nc.tol; ev1=[ev1 eplot(p)]; end 
end 
%% strategy 2, M_Vor, rlong=1, refine and retrig  
p=p0; p.sw.msw=0; p.sw.rlong=1; p=solfixpar(p); ev2=eplot(p); 
%
for i=1:nref;  p=refineX(p); p.X=setbcX(p); p=retrigX(p); %p=moveX(p,0.001,4);
 [p,res]=solfixpar(p); if res<p.nc.tol; ev2=[ev2 eplot(p)]; end 
end 