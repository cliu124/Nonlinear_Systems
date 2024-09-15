function qu=qjac2(p,u) 
X=p.X; M=getM(p,X); q1u=sum(M,1);
N=getN(p,p.X); if p.mpos; N=M*N; end 
qux=N(:,1)';  quy=N(:,2)'; 
qu=[q1u;qux;quy];