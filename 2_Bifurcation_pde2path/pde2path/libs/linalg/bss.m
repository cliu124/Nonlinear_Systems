function [x,y,p] = bss(A,b,c,d,f,g,p)
% BSS: bordered system solution after govaerts
% solve (A b; c' d)(x y)=(f g)
%  [x,y,p] = bss(A,b,c,d,f,g.p)
% 
% called by lssbss
% see also: lss, blss, lssbel, lssbelpo, lssbss, bel, belpo
[vs,p]=p.fuha.innerlss(A',c',p); dels=d-vs'*b;
[v,p]=p.fuha.innerlss(A,b,p); del=d-c*v;
y1=(g-vs'*f)/dels;
f1=f-b*y1; g1=g-d*y1; 
[w,p]=p.fuha.innerlss(A,f1,p); y2=(g1-c*w)/del;
x=w-v*y2; y=y1+y2;
end

