function p=spextend(p,spcont,prim)
% SPEXTEND: extend system with spectral part in case of spectral continuation
%  
%  p=spextend(p,spcont,prim)
%
% * p.sw.spcont is set to spcont in output p
% * p.nc.ilam in the output is augmented by prim at front
%
% See also spreduce, spcontini
    p.sw.spcont=spcont; p.nc.ilam=[prim;p.nc.ilam]; p.nu=2*p.nu; p.nc.neq=2*p.nc.neq; p.nc.nq=2*p.nc.nq+1; 
    if(p.sw.bcper~=0) % reset extended matrices in case of per. bcs
        p.mat.fill=blkdiag(p.mat.fill,p.mat.fill);
        p.mat.drop=blkdiag(p.mat.drop,p.mat.drop);
    end
end
