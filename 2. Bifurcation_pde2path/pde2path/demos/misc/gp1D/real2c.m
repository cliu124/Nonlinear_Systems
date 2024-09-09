function p=real2c(dir,pt,varargin) % convert real (scalar) GP to complex GP; 
% allow varargin to make changes in Newton method for initial steps: 
mod=0; if nargin>2; mod=1; newt=varargin{1}; tol=varargin{2}; end 
% load and plot soln to check that we have a decent (adapted) mesh 
p=loadp(dir,pt,[dir pt 'v']); plotsol(p,1,1,1,'pstyle','*'); 
M=p.mat.M; K=p.mat.K; par=p.u(p.nu+1:end); np=p.np;
p.mat.M=[M 0*M; 0*M M]; p.mat.K=[0*K K;-K 0*K]; % extend FEM matrices 
uv=[p.u(1:np); 0*p.u(1:np); par];  % extend solution 
p.nc.neq=2; p.nu=2*np; p.u=uv; p.sw.jac=0; p=resetc(p); p.nc.neig=20; 
p.fuha.sG=@vsG; p.fuha.outfu=@stanbra; p.com=1; % change rhs 
p.sw.verb=2; p.nc.amod=0; p.sw.bifcheck=0; % change some switches 
% and modify Newton method (if desired; -needed for some narrow peaks) 
if mod>0; oldtol=p.nc.tol; oldnewt=p.sw.newt; p.nc.tol=tol; p.sw.newt=newt; end 
p.sol.ds=1e-6; p=cont(p,2); % make two dummy steps to get on the branch 
if mod>0; p.nc.tol=oldtol; p.sw.newt=oldnewt; end % return to original Newton/tol 
plotspec(p,6);  title([dir '/' pt]); % plot the spectrum of cGP 
