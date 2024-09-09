function p=parabolinit(p,lx,ly,par,sw,nx) % parabol, just init mesh in 
% x-y-plane (w or wo symmetry controlled by sw), and map into R^3 
p=stanparam(p); p=stanparamX(p); p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; 
pde=stanpdeo2Db(lx,ly,nx,nx,sw); p.pdeo=pde; 
p.file.smod=1; p.sw.verb=2; p.nc.lammax=1e4; p.nc.lammin=-1e4; p.sw.bifcheck=2;  
po=pde.grid.p; x=po(1,:)'; y=po(2,:)'; 
a=par(1); b=par(2); z=x.^2/a^2+y.^2/b^2; % the 'exact solution' 
X=[x,y,z]; tri=pde.grid.t(1:3,:); tri=tri'; 
p.X=X; p.tri=tri; [p.tri,~]=orient_outward(p.X,p.tri);
p.np=size(p.X,1); p.nu=p.np; p.nt=size(p.tri,1); p.sol.xi=1/p.nu;
p.u=[zeros(p.np,1);par]; p.up=p.u; 
p.DIR=boundary_faces(p.tri); p.idx=unique(p.DIR); p.sw.nobdref=0; 
p.plot.auxdict={'H'}; p.tau=p.u(1:p.nu+1)'; p=setfn(p,'du'); 
p.nc.ilam=1; p.sw.para=0; p.sw.jac=0; p.sol.ds=1e-2; 
