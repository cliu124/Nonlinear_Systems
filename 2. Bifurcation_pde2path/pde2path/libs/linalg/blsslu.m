function [x,p]=blsslu(A,b,p)
% BLSSLU: check for LU in p.LU, and use that if up to date, otherwise update 
%
%  x=lsslu(A,au,p)
%
% See also blss, ilss
luok=1; havex=0; lucheck=2; d=0; 
try lsstol=p.nc.lsstol; catch; lsstol=1e-4; end 
try LU=p.mat.bLU; P=LU.P; Q=LU.Q; L=LU.L; U=LU.U; catch; luok=0; end 
if luok
   if size(L,1)~=size(A,1); luok=0; end % only check size 
   if luok; x=Q*(U\(L\(P*b))); havex=1; end 
   if lucheck==2 && havex;  % strict checking
     d=norm(A*x-b)/(1+norm(b)); if d>lsstol; luok=0; havex=0; end 
      %At=P'*(L*(U*Q')); d=max(max(abs(At-A))); 
   end
end
if ~luok; fprintf('blsslu: d=%g, new LU\n', d); 
    [L,U,P,Q]=lu(A); LU.L=L; LU.U=U; LU.P=P; LU.Q=Q; p.mat.bLU=LU; end 
if ~havex; x=Q*(U\(L\(P*b))); % d=norm(A*x-b)/(1+norm(b))
end 
end
   