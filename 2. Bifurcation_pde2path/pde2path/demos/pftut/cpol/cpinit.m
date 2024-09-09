function p=cpinit(p,lx,ly,nx,ny,par,ref,bpstyle,cut,h0) % AC on sphere coupled with diff-eqn in bulk
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; p.plot.cm='cool'; 
p.fuha.sG=@sGcp; p.fuha.sGjac=@sGcpjac; p.fuha.qf=@qfm; p.fuha.qfder=@qfmder; p.sw.jac=1; 
p.fuha.outfu=@cpbra; 
p.plot.auxdict={'R','k0','eps','gamma','s','m'}; p.plot.pstyle=-1; % userplot
pde1=spherepdeo(lx,ly,nx,ny,ref,1.2); % spherical coord. pdeo;  
p.pdeo=pde1; p.np=pde1.grid.nPoints; p.nu=p.np; 
p.u=0.3*ones(p.np,1); p.u=[p.u; par']; 
p=box2per(p,1); % pBC in x
p.nus=p.nu; p.nps=p.np; % store np and nu for sphere
u1=p.u(1:p.nus); uf=p.mat.fill*u1; 
figure(10); clf; R=par(1); p.pdeo.grid.spplot0(uf,R); % plot sphere 
po=getpte(p); xf=po(1,:)'; yf=po(2,:)'; 
x=p.mat.drop*xf; y=p.mat.drop*yf; 
i1=[]; i2=[]; ymin=min(y); ymax=max(y); 
for i=1:size(x)
    if y(i)==ymax; i1=[i1 i]; end % points adjacent to northpole 
    if y(i)==ymin; i2=[i2 i]; end % points adjacent to southpole 
end
p.i1=i1; p.i2=i2; p.nN=length(i1); p.nS=length(i2); 
pde2=sphere2ballpdeo(p,R,h0);
p.p2=pde2; p.npb=pde2.grid.nPoints; p.nub=p.npb;   % put ball-pde into p2 
% first p.nu pts in pde2 correpond to those in surface mesh, then N and S
p=oosetfemops(p); % now set up system matrices! 
ub=2.4*ones(p.npb,1); p.u=[u1; ub; par']; p.nc.ilam=2; p.usrlam=0:2:6; 
p.bpstyle=bpstyle; p.plot.cut=cut; p.plot.bpcmp=7;  
p.plot.sview=[10 20]; p.plot.cview=[15 40]; p.plot.pmod=2; 
p.nc.neig=10; p.nnc.eigref=-0.2; p.sw.bifcheck=2; p.nc.mu1=1; p.nc.mu2=0.01; 
%p=setilup(p,1e-4,200); p.sw.eigssol=0; p.sw.verb=2;  % not efficient! 
p.sol.ds=0.1; p.nc.dsmax=0.4; p.nc.ilam=6; % cont in mass m
p.S=0; % p.S=special sparsity-struc
p.sol.xi=1/(p.nu); 
