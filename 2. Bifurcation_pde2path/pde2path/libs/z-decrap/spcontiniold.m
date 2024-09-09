function p=spcontini(varargin)
% SPCONTINI: initialization for "spectral cont": continuation of zero
% eigenvalue for fold; MAY also work for branch-points (under symmetry) 
%
%  q=spcontini(dir,fname,inewpar)
%  q=spcontini(dir,fname,inewpar,newdir)    - set problem name to newdir
%  q=spcontini(dir,fname,inewpar,newdir,ds) - also set stepsize to ds
%
%  same with p, instead of dir/fname 
nnsw=0; resetsw=1; ds=0;
if ischar(varargin{1}) % old syntax
   dir=varargin{1}; fname=varargin{2}; inewpar=varargin{3}; 
   if nargin>=4; nname=varargin{4}; nnsw=1; end 
   if nargin==5; ds=varargin{5}; end 
   p=loadp(dir,fname);    %generates also p.mat.drop and p.mat.fill
else p=varargin{1}; inewpar=varargin{2}; 
   if nargin>=3; nname=varargin{3}; if ~isempty(nname); nnsw=1; end; end 
   if nargin==4; ds=varargin{4}; end 
end 
if (p.sol.ptype~=2 && p.sol.ptype~=1)
    fprintf('\nWarning: not a fold point. Set p.sw.spcont to 2 (fold point) by hand;\n'); 
    fprintf('the Newton loop may converge to the desired point\n'); 
end
fprintf('Initializing ... \n'); 
if (p.sw.jac~=1) fprintf('\nWarning: Computation may be slow as all derivatives are computed numerically.\n'); end
if (p.sw.qjac~=1) fprintf('\nWarning: Computation may be slow as some pde derivatives are computed numerically.\n'); end
% Compute leading eigenvector
M=getM(p); M=[[M zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) eye(p.nc.nq)]]; Gu=getGu(p,p.u,resi(p,p.u));
%double the size of p.mat.drop and p.mat.fill  
r=size(p.mat.drop,1); s=size(p.mat.drop,2);
if r>1;  p.mat.drop = [[p.mat.drop zeros(r,s)]; [zeros(r,s) p.mat.drop]]; 
  p.mat.fill = [[p.mat.fill zeros(s,r)]; [zeros(s,r) p.mat.fill]]; end
evopts.disp=0; [phi1,mu1]=eigs(Gu,M,1,0,evopts);
fprintf('Approx. zero eigenvalue=%g.\n',mu1); phi1=phi1/norm(phi1); % normalize eig-fct. 

% Create new initial vector: format is u=[updeold; updelin; uauxold; auxlin]
% Here auxlin has the length of active aux. vars. except primary par.
p.naux=length(p.u)-p.nu; % store number of auxiliary variables 
p.u=[p.u(1:p.nu); phi1(1:p.nu); p.u(p.nu+1:p.nu+p.naux); phi1(p.nu+1:end)];
% set new active variables including those for linear part (length=p.nc.nq+2+p.nc.nq)
p.nc.ilam=[inewpar; reshape(p.nc.ilam,length(p.nc.ilam),1); [p.naux+1:p.naux+p.nc.nq]'];

% set new equation numbers
p.nc.neq=2*p.nc.neq; p.nu=2*p.nu; p.nc.nq=2*p.nc.nq+1; % HU 
p.sw.spcont=2; p.sol.restart=1; 
fprintf('New active parameters and their values:\n'); printaux(p,1)
% other data initialization
if ds~=0; p.sol.ds=ds; end
if resetsw; p=resetc(p); end % reset counters and branch
%if(~isempty(varargin)); [p,ok]=setfn(p,varargin{1}); if ok~=1; q=p; return; end
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
else fprintf('warning: file name prefixes unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1; p.sw.spcalc=0;
fname=[p.file.pname,sprintf('%i',p.file.count),'.mat'];save(fname,'p'); 