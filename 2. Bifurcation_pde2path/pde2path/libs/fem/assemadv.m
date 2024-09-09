function K=assemadv(p,t,b)
% ASSEMADV: assemble K for advective terms. 
%
% K=assemadv(p,t,b)
% b=b(i,k), i=1..2N^2, k=1,..,nt is row-vector on triangles 
% p=mesh-points, t=triangles, e=edges, N=# of eqns/unknowns (dim of system)
%
% See also assemadv1, or (OOPDE) fem.convection
N=sqrt(size(b,1)/2); 
m1=1;m2=2; np=size(p,2); ks=sparse(N*np,N*np,0); 
for l=1:N % column blocks in K 
    for k=1:N % row-blocks in K
        if(any(b(m1:m2,:))) ks1=assemadv1(p,t,b(m1:m2,:)); 
          [ii,jj,kss]=find(ks1);
          ks=ks+sparse(ii+(k-1)*np,jj+(l-1)*np,kss,N*np,N*np);
        end
        m1=m1+2; m2=m2+2; 
    end
end
K=ks;
