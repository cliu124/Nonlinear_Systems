function p=loadp(dir,fname,varargin)
% load point from dir/fname.mat
noa=nargin-2; % number of opt arguments
ffname=[dir '/' fname '.mat'];
s=load(ffname,'p'); p=s.p; fprintf('lam=%g\n',getlam(p));
if noa>0; p=setfn(p,varargin{1}); end;

% set operators that were not saved in p
if(0)
    if(p.sw.bcper==3)
        [p.mat.fill, p.mat.drop, p.nu]=getTorOp(p);
    elseif(p.sw.bcper==1 || p.sw.bcper==2)
        [p.mat.fill, p.mat.drop, p.nu]=getCylOp(p);
    end
end
[p.mat.fill, p.mat.drop, p.nu]=getPerOp(p);

if(p.sw.sfem~=0 || p.sw.spcalc~=0) p=setfemops(p); end
if(p.pffac<0)
    onu=p.nu; par=p.u(p.nu+1:end); p=dropp(p);
    lu=length(p.tau) ;
    p.tau=[p.tau(1:lu/2+p.pfn-1); p.tau(lu/2+p.pfn+1:lu)];
end
p.mat.pot=pot(p);p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot);
