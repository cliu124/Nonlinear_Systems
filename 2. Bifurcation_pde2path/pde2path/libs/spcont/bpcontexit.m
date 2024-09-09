function p=bpcontexit(varargin)
% BPCONTEXIT: exit BP continuation
%
%  p=bpcontexit(dir,fname)
%  p=bpcontexit(dir,fname,newdir) - set problem name to newdir
%  p=bpcontexit(p,newdir) also allowed 
nnsw=0; 
if ischar(varargin{1}); dir=varargin{1}; fname=varargin{2}; 
   if nargin>2; nname=varargin{3}; nnsw=1; end % nnsw='new-name-switch' 
   ffname=[dir '/' fname]; s=load(ffname,'p'); p=s.p;
   if nargin>3; p.fuha.outfu=varargin{4}; end 
else p=varargin{1}; if nargin==2; nname=varargin{2}; nnsw=1; end 
end 
if(p.sw.spcont==0) fprintf('Not a fold or branch point continuation problem.\n'); q=p; return; end
p.nu=p.nu/2; p.nc.neq=p.nc.neq/2; p.nc.nq=(p.nc.nq-2)/2; % set regular case sizes
p.sol.ptype=p.sw.spcont; p.sw.spcont=0;
p.u=[p.u(1:p.nu);p.u(2*p.nu+1:2*p.nu+p.naux)]; % non-fold/branch point variables
p.nc.ilam=p.nc.ilam(2:2+p.nc.nq); % remove current primary parameter
p.tau=[p.tau(1:p.nu);p.tau(2*p.nu+1:2*p.nu+p.nc.nq);0];
fprintf('Exiting branch point continuation. Now\n'); printaux(p,1); 
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); % generate p.mat.fill, p.mat.drop (not saved in data files)
if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end
% other data initialization
if nnsw==1; p.file.bcount=1; p.file.count=0; end 
brout=p.fuha.outfu(p,p.u); % userfu to append to bif-branches
try p.branch=[p.branch [bradat(p); brout]]; % put on branch
catch p.branch=[bradat(p); brout];
end 
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
   p.branch=p.branch(1:length(bradat(p))+length(p.fuha.outfu(p,p.u)),end); % cut off branch to new length 
else fprintf('warning: problem directory unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1;
p.fuha.savefu(p); p.file.count=p.file.count+1; 
% cont in both directions to have tangent for swibra
p.sol.ds=-1e-6; % p.plot.bpcmp=0; 
p.sw.spcalc=0; p.nc.tol=1e-8; bifsw=p.sw.bifcheck; p.sw.bifcheck=0; 
p=cont(p,2); p.sol.ds=-p.sol.ds; p=resetc(p); p=cont(p,1); p.sw.bifcheck=bifsw;  