function p=enninit(nx,par) % Enneper init 
p=stanparam(); p=stanparamX(p); p.sw.bifcheck=1; % set stanparam, overwrite some
p.fuha.sG=@sGenn; p.fuha.outfu=@cmcbra; p.fuha.e2rs=@e2rsbdry; p.sw.jac=0;
al=par(4); pde=diskpdeo2(al,nx,round(nx/2)); % disk preimage discretization, pde-object 
% not stored, only p.DBC, p.tri and p.X (generated below) used subsequently
p.np=pde.grid.nPoints; p.nu=p.np; p.nt=pde.grid.nElements;  % store dimensions
po=pde.grid.p; x=po(1,:)'; y=po(2,:)'; p.sol.xi=1/p.nu; p.nc.neq=1; p.sw.sfem=-1; 
u=0*ones(p.np,1); p.u=[u; par]; r=sqrt(x.^2+y.^2); th=angle(x+1i*y); % set ICs 
p.X=[r.*cos(th)-r.^3.*cos(3*th)/3, -r.*sin(th)-r.^3.*sin(3*th)/3, r.^2.*cos(2*th)];
p.tri=pde.grid.t(1:3,:)'; % store initial X and tri 
p.DIR=pde.grid.e(1:2,:)'; p.idx=unique(p.DIR(:)); % edges, and points for BCs 
p=oosetfemops(p); p.plot.auxdict={'H','V','A','alpha','k'}; p.plot.pstyle=-1; % userplot on 
p.u(p.nu+2)=getV(p,p.u); p.u(p.nu+3)=getA(p,p.u); % get initial vol & area
p.sw.verb=2; p.plot.bpcmp=7; p.file.smod=1; 
p.nc.lammin=0; p.nc.lammax=2; p.th=th(p.idx); 
p.sw.rlong=1; p.sw.nobdref=0; p.fuha.e2rs=@e2rsbdry; 