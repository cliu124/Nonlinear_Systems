function [u,res,iter,Gu,Glam]=fsol(p,u)
% FSOLEXT: replaces nloopext for using fsol.
%
%  [u,res,iter,Gu,Glam]=fsol(p,u)
% Important settings: p.fsol.opt, p.fsol.tol, p.fsol.imax, p.fsol.meth, p.fsol.disp
%
% See also nloop, fsolext, fsolextf, stanparam
global pfs;
pfs=p; ua=u2au(p,u); 
opt=optimset(p.fsol.opt,'TolFun',p.fsol.tol,'MaxIter',p.fsol.imax);
if p.fsol.meth==1; opt=optimset(opt,'Algorithm','trust-region-reflective');
else opt=optimset(opt,'Algorithm','Levenberg-Marquardt'); end
switch p.fsol.disp
    case 0; opt=optimset(opt,'Display','off');
    case 1; opt=optimset(opt,'Display','final');
    case 2; opt=optimset(opt,'Display','notify');
    case 3; opt=optimset(opt,'Display','iter');
end
[x,r,exitflag,output] = fsolve(@fsolvef,ua,opt);
u=au2u(p,x); res=norm(r,p.sw.norm); 
iter=output.iterations; [Gu,Glam]=getder(p,u,r);
end
