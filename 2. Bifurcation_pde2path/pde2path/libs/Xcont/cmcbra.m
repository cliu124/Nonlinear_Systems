function out=cmcbra(p,u) 
% cmcbra: template output for geometric problems
%
% out=cmcbra(p,u); 
% 
% after pars we put V, A, E, h1, h2, A1, A2, q1, q2,  
% h1=max(edgelengths), h2=min(edgelengths), 
% A1=max(triangle areas), A2=min(triangle areas) 
% q1=max(h/r), q2=max(R/r),   r=incircle, R=outcircle, 
par=u(p.nu+1:end); X=p.X; V=getV(p,u); Ai=getA(p,u); E=Ai+par(1)*V; 
mqd=meshqdat(p); % max(h/r), max(A), min(A), max(h), min(h) 
out=[par; V;Ai;E; mqd]; 