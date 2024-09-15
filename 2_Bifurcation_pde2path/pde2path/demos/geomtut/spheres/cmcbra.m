function out=cmcbra(p,u) 
% mod of cmcbra; after pars we put V, A, E, ***r***, mqd 
par=u(p.nu+1:end); X=p.X; V=getV(p,u); Ai=getA(p,u); E=Ai+par(1)*V; 
n=sqrt(dot(X,X,2)); r=max(n);  
mqd=meshqdat(p); % max(h/r), max(A), min(A), max(h), min(h) 
out=[par; V;Ai;E; r; mqd]; 