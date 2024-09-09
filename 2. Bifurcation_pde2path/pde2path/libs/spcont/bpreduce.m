function [p,spcont,prim]=bpreduce(p)
% BPREDUCE: reduce to original system in case of BPcontinuation
spcont=p.sw.spcont; p.sw.spcont=0; prim=p.nc.ilam(1); 
p.nu=p.nu/2; p.nc.neq=p.nc.neq/2; p.nc.nq=(p.nc.nq-2)/2; % set regular case sizes
p.nc.ilam=p.nc.ilam(2:end-1-p.nc.nq);
if(p.sw.bcper~=0) % submatrices in case of periodic bcs
     p.mat.fill=p.mat.fill(1:p.nc.neq*p.np,1:p.nu);
     p.mat.drop=p.mat.drop(1:p.nu,1:p.nc.neq*p.np);
end