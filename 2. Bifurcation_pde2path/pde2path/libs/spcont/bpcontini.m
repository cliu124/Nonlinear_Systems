function p=bpcontini(varargin)
% BPCONTINI: initialization for BP continuation 
%
% p=bpcontini(dir,fname,newpar) 
% p=bpcontini(dir,fname,newpar,newdir,ds) 
% 
% newpar is parameter for BPcont, old-parameter(s) automatically appended 
%  same with p, instead of dir/fname 
nnsw=0; resetsw=1; ds=0;
if ischar(varargin{1}) % old syntax
   dir=varargin{1}; fname=varargin{2}; inewpar=varargin{3}; 
   if nargin>=4; nname=varargin{4}; nnsw=1; end; if nargin==5; ds=varargin{5}; end 
   p=loadp(dir,fname);    %generates also p.mat.drop and p.mat.fill
else p=varargin{1}; inewpar=varargin{2}; 
   if nargin>=3; nname=varargin{3}; if ~isempty(nname); nnsw=1; end; end 
   if nargin==4; ds=varargin{4}; end 
end 
if p.sol.ptype~=1 
    fprintf('\nWarning: not a branch point. Set p.sw.spcont to 1 (branch point) by hand.\n'); 
    fprintf('The Newton loop might converge to the desired point\n'); 
    return; 
end
fprintf('Initializing branch point continuation.\n'); 
if p.nc.nq>0; p.fuha.quupsi=@quupsi; end % set default function handle for q_uu
if (p.sw.jac~=1); 
    fprintf('\nWarning: Computation may be slow as PDE derivatives are computed numerically.\n'); end
if (p.nc.nq>0 && p.sw.qjac~=1); 
    fprintf('\nWarning: Computation may be slow as constraint derivatives are computed numerically.\n'); end
% Compute 0 eigenvector of Gu^T
M=getM(p); M=[[M zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) 0*speye(p.nc.nq)]]; 
Gu=getGu(p,p.u,resi(p,p.u));
evopts.disp=0; [psi1,mu1]=eigs(Gu',M,1,0,evopts); 
fprintf('Approx. zero eigenvalue=%g.\n',mu1); psi1=psi1/norm(psi1); % normalize EVec
% Create new initial vector:  u=[upde; psi; uauxold; auxlin; mu]
% Here auxlin has the length of active auxiliary vars. except primary par.
p.naux=length(p.u)-p.nu; % store number of auxiliary variables 
p.u=[p.u(1:p.nu); psi1(1:p.nu); p.u(p.nu+1:p.nu+p.naux); mu1; psi1(p.nu+1:end)]; 
% set new active variables including those for linear part (length=p.nc.nq+2+p.nc.nq)
p.nc.ilam=[inewpar,p.nc.ilam,p.naux+1:p.naux+1+p.nc.nq];
%           new        old     mu and psi_q (aux-compos of adj.eigenvector)                
p=doubledrop(p); % double the size of p.mat.drop and p.mat.fill
% set new equation numbers
p.nc.neq=2*p.nc.neq; p.nu=2*p.nu; p.nc.nq=2*p.nc.nq+2; 
p.sw.spcont=1; p.sol.restart=1;
fprintf('New active parameters and their values:\n'); printaux(p,1)
% other data initialization
if ds~=0; p.sol.ds=ds; end
if resetsw; p=resetc(p); end % reset counters and branch
%if(~isempty(varargin)); [p,ok]=setfn(p,varargin{1}); if ok~=1; q=p; return; end
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
else fprintf('warning: file name prefixes unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1; p.sw.spcalc=0;
p.plot.bpcmp=p.nc.ilam(2); % old primary! 
p.sw.spjac=1; p.fuha.spjac=@bpjac; 
p.sw.bifcheck=0; p.sw.spcalc=0; p.nc.del=1e-4; % set some switches 
fname=[p.file.pname,sprintf('%i',p.file.count),'.mat'];save(fname,'p'); 