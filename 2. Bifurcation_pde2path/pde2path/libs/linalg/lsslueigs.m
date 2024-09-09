function [x,p]=lsslueigs(A,b,p)
% LSSLUeigs: check for LU in p.LU, and use that if up to date
global p2pglob; 
luok=1; havex=0; lucheck=2; d=0; 
try lsstol=p.nc.lsstol; catch; lsstol=1e-6; end 
try; LU=p2pglob.LU; P=LU.P; Q=LU.Q; L=LU.L; U=LU.U; catch; luok=0; end 
if luok
   if size(L,1)~=size(A,1); luok=0; end % only check size 
   if luok==1 && lucheck==1
       At=P'*(L*(U*Q')); d=max(max(abs(At-A))); 
       if d>lsstol; luok=0; end 
   end
   if luok; x=Q*(U\(L\(P*b))); havex=1; end 
   if lucheck==2 && havex  % check residual
      d=norm(A*x-b)/(1+norm(b));  
      if d>lsstol; luok=0; havex=0; end      
   end
end
if ~luok; if p.sw.verb>1; fprintf('lsslueigs: d=%g, new LU\n', full(d)); end 
    if ~issparse(A); A=sparse(A); end 
    [L,U,P,Q]=lu(A); LU.L=L; LU.U=U; LU.P=P; LU.Q=Q; p2pglob.LU=LU; 
end 
if ~havex; x=Q*(U\(L\(P*b)));  
end 
end 