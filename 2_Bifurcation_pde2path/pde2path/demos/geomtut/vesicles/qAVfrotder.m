function qu=qAVfrotder(p,u) % u-derivative of A,V,full transl. and rot. PCs 
par=u(p.nu+1:end); H=u(p.np+1:p.nu); u=u(1:p.np); X0=p.X; N=getN(p,X0); 
%Rotation
grX1dc=cross(p.om.*ones(p.np,3),X0,2); uph1=dot(grX1dc,N,2);
grX1dc=cross(p.rh.*ones(p.np,3),X0,2); uph2=dot(grX1dc,N,2);
grX1dc=cross(p.w.*ones(p.np,3),X0,2); uph3=dot(grX1dc,N,2);
%Translation
dsx=N(:,1); dsy=N(:,2); dsz=N(:,3); 
M=getM(p,X0); M=M(1:p.np,1:p.np); 
%M1=getM(p,X0+u.*N); M1=M1(1:p.np,1:p.np);
M1=massmatrix(X0+u.*N,p.tri,'voronoi'); 
qu=[(-2*M1*H)' 0*dsx'; 
    -sum(M1,1) 0*dsx'; 
    (M*dsx)' 0*dsx'; 
    (M*dsy)' 0*dsx';  
    (M*dsz)' 0*dsx';...
    (M*uph1)' 0*dsx'; 
    (M*uph2)' 0*dsx'; 
    (M*uph3)' 0*dsx';];
