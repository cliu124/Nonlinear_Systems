function [f,fail]=nleq1link(u,FLAG,p) 
% nleq1link: wrapper for resi and Jac for NLEQ1 Newton loops
fail=0; %par=p.u(p.nu+1:end); up=[u; par]; 
switch FLAG
    case ''; f=nleq1resi(p,u); % size(f) %f=resi(p,up); 
    case 'jacobian'; if p.sw.jac==0; r=resi(p,up); else; r=0; end 
        f=nleq1getGu(p,u,r);  % size(f), p.nu, %f=getGu(p,up,r);
end
