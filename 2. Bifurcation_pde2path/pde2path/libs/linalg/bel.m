function [x,y,p]=bel(A,b,c,d,f,g,p)
% BEL: bordered elimination, called by lssbel 
%
%  [x,y,p]=bel(A,b,c,d,f,g,p)
% 
% called by belpi
% see also: lss, blss, lssbel, blssbel, belpi
[v,p]=p.fuha.innerlss(A,b,p); 
if norm(f,'inf')>1e-10; [w,p]=p.fuha.innerlss(A,f,p); 
else w=zeros(size(f,1),1); end 
del=d-c*v; y=del\(g-c*w); x=w-v*y;

