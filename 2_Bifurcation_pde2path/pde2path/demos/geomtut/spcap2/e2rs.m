function [p,idx]=e2rs(p,u)  % elements2refine via triangle size on X 
X=p.X; Xx=p.mat.Dx*X; Xy=p.mat.Dy*X; n=cross(Xx,Xy,2);
E=dot(n,n,2); E=p.mat.p2c*E; % map to elements 
[E,idx]=sort(E,'descend'); idx=idx(1:round(p.nc.sig*length(E))); 
p.sol.err=max(max(E)); 