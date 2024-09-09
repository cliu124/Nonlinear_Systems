function [di,d2]=tadev(p)
% tadev: return target deviation in sup and l2 norm 
uend=p.cp.u(:,end); nu=length(uend); uh1=p.oc.u1(1:nu); di=max(abs(uh1-uend)); 
d2=sqrt(norm(uh1-uend,2)^2/nu);
fprintf('target sup-deviation=%g, target l2-deviation=%g\n',di,d2); 