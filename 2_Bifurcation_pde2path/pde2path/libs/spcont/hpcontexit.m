function p=hpcontexit(varargin)
% HPCONTEXIT: exit Hopf point continuation 
%
%  p=hpcontexit(dir,fname)
%  p=hpcontexit(dir,fname,newdir)    - set problem name to newdir
%  p=hpcontexit(dir,fname,newdir,outfu) - set p.fuha.outfu=outfu; 
nnsw=0; % new name switch 
if ischar(varargin{1}); dir=varargin{1}; fname=varargin{2}; 
   if nargin>2; nname=varargin{3}; nnsw=1; end 
   ffname=[dir '/' fname]; s=load(ffname,'p'); p=s.p;
   if nargin>3; p.fuha.outfu=varargin{4}; end 
else p=varargin{1}; if nargin==2; nname=varargin{2}; nnsw=1; end 
end 
if(p.sw.spcont~=3) fprintf('Not a HP continuation problem.\n'); return; end
p.nu=p.nu/3; p.nc.neq=p.nc.neq/3; p.nc.nq=(p.nc.nq-2)/3; % set regular case sizes
p.sol.ptype=p.sw.spcont; p.sw.spcont=0;
p.u=[p.u(1:p.nu);p.u(3*p.nu+1:3*p.nu+p.naux)]; % non-fold/branch point variables
p.nc.ilam=p.nc.ilam(2:2+p.nc.nq); % remove current primary parameter
p.tau=[p.tau(1:p.nu);p.tau(3*p.nu+1:3*p.nu+1+p.nc.nq);0];
fprintf('Exiting HP continuation. Now\n'); printaux(p,1); 
%generate p.mat.fill, p.mat.drop (these are not saved in data files)
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p); 
if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end
% other data initialization
if nnsw==1; p.file.bcount=1; p.file.count=0; end 
brout=p.fuha.outfu(p,p.u); % userfu to append to bif-branches
try; p.branch=[p.branch [bradat(p); brout]]; % put on branch
catch; p.branch=[bradat(p); brout];
end 
if nnsw; [p,ok]=setfn(p,nname); if ok~=1; q=p; return; end
   p.branch=p.branch(1:length(bradat(p))+length(p.fuha.outfu(p,p.u)),end); % cut off branch to new length 
else fprintf('warning: problem directory unchanged.\n'); end
p.sol.deta=0; p.sol.ineg=-1;
p.fuha.savefu(p); p.file.count=p.file.count+1; 
%fname=[p.file.pname,sprintf('%i',p.file.count),'.mat']; save(fname,'p');
%p.file.count=2;p.file.bcount=1;