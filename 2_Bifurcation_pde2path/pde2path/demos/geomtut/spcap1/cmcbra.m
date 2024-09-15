function out=cmcbra(p,u)  % mod of template cmcbra with error-output 
par=u(p.nu+1:end); X=p.X; V=getV(p,u); A=getA(p,u); E=A+par(1)*V; 
Z2=max(X(:,3)); mqdat=meshqdat(p); 
% augment lib-cmcbra by specialized output, here: errors  
svp=sqrt(pi^2+9*V^2); % shorthand 
HV=-pi^(1/3)*(3*V+svp-pi^(2/3)*(svp-3*V)^(1/3))/(svp*(3*V+svp)^(1/3)); % analytical H 
Z2V=((3*V+svp)^(2/3)-pi^(2/3))/(pi^(1/3)*(3*V+svp)^(1/3)); % analytical height  
M=getM(p,X); LB=cotmatrix(X,p.tri); N1=getN(p,X); H1=0.5*dot(LB*X,N1,2); 
q=boundary_faces(p.tri); ids=unique([q(:,1);q(:,2)]); idb=setdiff([1:p.np],ids); % bulk indices 
H2=M\H1; eHinf=max(abs(HV-H2(idb)))/abs(HV); % max relative error in H (in bulk) 
eH2=full(sqrt(sum(M(idb:idb)*(HV-H2(idb)).^2)))/abs(HV); % L2-error (over bulk) 
out=[par;  V;A;E; mqdat; Z2; eH2; eHinf; abs(Z2V-Z2)]; 
%   #par  +1      4-8    9   10   11       12 