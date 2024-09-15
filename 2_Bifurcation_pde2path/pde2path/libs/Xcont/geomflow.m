function [X,t,ts]=geomflow(p,t,ts,dt,ns,nplot) 
% geomflow: explicit Euler stepping X_{n+1}=X_n+dt*f(X_n)*N, 
%
% f interfaced in p.fuha.flowf
par=p.u(p.nu+1:end); X=p.X; 
for i=0:ns-1; 
    f=p.fuha.flowf(p,X); % the (scalar) rhs; 
    if mod(i,nplot)==0;  % check for plotting and time-series output
       p.up=[f;par]; p.X=X; p.t=t; pplot(p); 
       A=getA(p,0*p.u); V=getV(p,0*p.u); mq=meshqdat(p); 
       ts=[ts [t; A; V;mq]]; fprintf(' t=%g, V=%g,  ',t,V); 
    end 
    N=getN(p,X); X=X+dt*f.*N; t=t+dt; % flow 
end 
p.up=[f;par]; p.X=X; p.t=t; pplot(p); 
A=getA(p,0*p.u); V=getV(p,0*p.u); mq=meshqdat(p); 
ts=[ts [t; A; V;mq]]; fprintf(' t=%g, V=%g,  ',t,V); 
fprintf('\n'); 