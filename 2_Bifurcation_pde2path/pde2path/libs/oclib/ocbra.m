function out=ocbra(p,u) 
% ocbra: standard output to bifurcation diagram for OC 
%
% out=ocbra(p,u) 
[po,tr,e]=getpte(p); np=p.np; M=p.mat.M(1:np,1:np); 
out=[max(u(1:p.np)); min(u(1:p.np)); 
    sqrt(sum(M*(u(1:np).^2),1)/p.vol); jca(p,u)];
end