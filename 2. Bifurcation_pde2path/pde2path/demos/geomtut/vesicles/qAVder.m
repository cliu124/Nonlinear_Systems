function qu=qAVder(p,u) % A,V and 3 translational PCs 
H=u(p.np+1:p.nu); u=u(1:p.np); X0=p.X; N=getN(p,X0); 
dsx=N(:,1); dsy=N(:,2); dsz=N(:,3); 
M=getM(p,X0); M=M(1:p.np,1:p.np); 
M1=massmatrix(X0+u.*N,p.tri,'full'); % seems better for vol constraint
qu=[(-2*M*H)' 0*dsx'; 
    -sum(M1,1) 0*dsx'; 
    (M*dsx)' 0*dsx'; 
    (M*dsy)' 0*dsx';  
    (M*dsz)' 0*dsx';];

