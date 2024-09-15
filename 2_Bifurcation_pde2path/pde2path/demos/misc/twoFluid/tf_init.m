function p=tf_init(p,nx,ny,lx,ly) 
% standard settings, file names and function handles
p=stanparam(p); dir=sprintf('%s',inputname(1)); p=setfn(p,dir); screenlayout(p);
% PDE
p.fuha.G=@tf_G; p.fuha.Gjac=@tf_Gjac; p.nc.neq=3;   

% domain
p.mesh.geo=rec(lx/2.,ly/2.); 

% BC: Dirichlet left/right and Neumann top/bottom  
qd = 10^3*[[1;0;0] [0;1;0] [0;0;1]]; % Dirichlet 
qn = zeros(p.nc.neq); gv = zeros(p.nc.neq,1);
bc = gnbc(p.nc.neq,qn,gv,qd,gv,qn,gv,qd,gv);
p.fuha.bc=@(p,u) bc; p.fuha.bcjac=@(p,u) bc;
p=stanmesh(p,nx,ny); p=setbmesh(p); 

% continuation settings
p.nc.ilam=1;
p.sol.ds=-0.001; p.nc.dsmax=0.1; p.nc.lammax=2; p.nc.lammin=0.1;
p.sol.xi=1/p.np; p.nc.lamdtol=0.8; 

% initial point and parameters
p.u=0*ones(p.nu,1); % NB: for cylinder bc u is the reduced vector, while the mesh is for full
delT=0.159;%(2*L2^2/L1)/(4+L2^2/L1^2)/pi^2; 
nu=9e-4; om=delT*pi/ly; s=-om*ly/(2*pi); 
par(1)=delT; par(2)=nu; par(3)=lx; par(4)=ly; par(5)=s;  
p.u = [p.u; par'];

% switch to cylinder bc with adaption of initial guess
p=box2per(p,2);
