function f=mcff(p,X)
% mcff: rhs f (scalar) for MCF
N=getN(p,X); M=getM(p,X); LB=cotmatrix(X,p.tri); 
f=0.5*dot(LB*X,N,2); f=M\f; f(p.idx)=0; 
