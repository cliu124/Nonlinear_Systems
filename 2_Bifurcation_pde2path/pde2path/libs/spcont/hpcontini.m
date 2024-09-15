function p=hpcontini(varargin)
% HPCONTINI: initialization for Hopf-point continuation 
%  p=hpcontini(dir,fname,inewpar)
%  p=hpcontini(dir,fname,inewpar,newdir)    - set problem name to newdir
%  p=hpcontini(dir,fname,inewpar,newdir,ds) - also set stepsize to ds
%  p=hpcontini(dir,fname,inewpar,newdir,ds,aux) - add.arguments in aux
%
%  same with p, instead of dir/fname 
nnsw=0; resetsw=1; ds=0; 
if ischar(varargin{1}) % old syntax
   dir=varargin{1}; fname=varargin{2}; inewpar=varargin{3}; 
   if nargin>3; nname=varargin{4}; nnsw=1; end 
   if nargin>4; ds=varargin{5}; end 
   if nargin>5; aux=varargin{6}; end 
   p=loadp(dir,fname);    %generates also p.mat.drop and p.mat.fill
else p=varargin{1}; inewpar=varargin{2}; 
   if nargin>2; nname=varargin{3}; if ~isempty(nname); nnsw=1; end; end 
   if nargin>3; ds=varargin{4}; end 
   if nargin>4; aux=varargin{5}; end 
end 
if (p.sol.ptype~=3 && p.sol.ptype~=1)
    fprintf('\nWarning: not a HP. Set p.sw.spcont to 3 (HP) by hand;\n'); 
    fprintf('the Newton loop may converge to the desired point\n'); 
end
if p.nc.nq>0; p.fuha.quuphir=@quuphir; p.fuha.quuphii=@quuphii; end % set default function handle for q_uu
fprintf('Initializing ... \n'); 
% Compute leading eigenvector
M=getM(p); M=[[M zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) 0*speye(p.nc.nq)]]; 
Gu=getGu(p,p.u,resi(p,p.u));
%om0=p.nc.eigref(2); p.sol.muv, 
try eignr=aux.eignr; catch; eignr=length(p.nc.eigref); end, eignr
om0=p.sol.muv(eignr,1); % om for finding Evec
evopts.disp=0; [phi1,mu1]=eigs(Gu,M,1,om0,evopts); 
fprintf('real(mu)=%g, imag(mu)=%g.\n',real(mu1),imag(mu1)); phi1=phi1/norm(phi1); % normalize eig-fct. 
phir=real(phi1); phii=imag(phi1); om=imag(mu1); p.c=phir'./norm(phir,2).^2; 
p.naux=length(p.u)-p.nu; % store (old) number of auxiliary variables 
% Create new initial vector: u=[u; phi_r; phi_i; pars; om; phi_rq; phi_iq]
p.u=[p.u(1:p.nu); phir(1:p.nu); phii(1:p.nu); p.u(p.nu+1:p.nu+p.naux); 
     om; phir(p.nu+1:end); phii(p.nu+1:end)];
% set new active variables including those for linear part (length=p.nc.nq+2+p.nc.nq)
p.nc.ilam=[inewpar,reshape(p.nc.ilam,1,length(p.nc.ilam)),p.naux+1:p.naux+1+2*p.nc.nq];
% set new equation numbers
p.nc.neq=3*p.nc.neq; p.nu=3*p.nu; p.nc.nq=3*p.nc.nq+2; 
p.sw.spcont=3; p.sol.restart=1;
p=tripledrop(p); % triple the size of p.mat.drop and p.mat.fill
fprintf('New active parameters and their values:\n'); printaux(p,1)
if ds~=0; p.sol.ds=ds; end % other data initialization
if resetsw; p=resetc(p); end % reset counters and branch
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
else fprintf('warning: file name prefixes unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1; p.sw.spcalc=0;
fname=[p.file.pname,sprintf('%i',p.file.count),'.mat'];save(fname,'p'); 
p.sw.spjac=1; p.fuha.spjac=@hpjac; 
p.sw.bifcheck=0; p.sw.spcalc=0; p.nc.del=1e-4; % set some switches 
