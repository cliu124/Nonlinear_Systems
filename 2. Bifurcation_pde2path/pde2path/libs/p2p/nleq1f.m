function [f,fail]=nleq1f(au,FLAG,p) 
% nleq1f: wrapper for resi and Jac for NLEQ1 Newton loops
fail=0; u=au2u(p,au); 
switch FLAG
    case ''; f=nleq1resi(p,u); % size(f) %f=resi(p,up); 
    case 'jacobian'; % if p.sw.jac==0; r=resi(p,u); else; r=0; end 
        r=resi(p,u); 
        f=nleq1getGu(p,u,r);  % size(f), p.nu, %f=getGu(p,up,r);
end
