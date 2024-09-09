function f=vpmcff(p,X) 
% vpmcff: rhs f for volume preserving MCF, explicit, and BCs ignored; 
% reasonable conservation of V, and convergence to steady state 
N=getN(p,X); M=getM(p,X); LB=cotmatrix(X,p.tri); H=0.5*dot(LB*X,N,2); 
H1=M\H; Hb=sum(H,1); % H=weak, H1=nodal values, Hb=integral 
A=doublearea(X,p.tri)/2; A=sum(A);  Hb=Hb/A; f=H1-Hb;