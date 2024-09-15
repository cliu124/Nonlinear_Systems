function [p,idx]=mye2rs(p,u) % ad-hoc triangle selcetion by usr, here: 
% select p.nr triangles with points closest to p.xrp; 
p.sol.err=0; [po,t,e]=getpte(p); 
[ps,idx]=sort(abs(po-p.xrp),'ascend'); pnr=idx(1:p.ntr); 
[lia,idx]=ismember(pnr,t(1,:));

