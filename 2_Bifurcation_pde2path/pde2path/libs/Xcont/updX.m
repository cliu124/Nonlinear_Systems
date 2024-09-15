function  [p,u]=updX(p,u)  
% updX: template for update p.X, return u in p.up, set u to zero
%
% [p,u]=updX(p,u)  
N=getN(p,p.X); np=p.np; if p.sw.Xcont<2; p.up=u; end
p.X=p.X+(p.mat.fill*u(1:np)).*N; u(1:np)=0; 