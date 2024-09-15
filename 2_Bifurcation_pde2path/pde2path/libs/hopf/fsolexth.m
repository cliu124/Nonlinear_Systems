function [y,T,lam,res,iter,A,p,jac]=fsolext(p,y1,T1,lam1,ds)
% fsolext: interface to call fsolve for Hopf in arclength parametrization
% 
global pfs;
pfs=p; 
nu=p.nu; tl=p.hopf.tl; na=nu;  yv=reshape(y1(1:na,:),na*tl,1); % make y-vector
ua=[yv; T1; lam1; p.u(p.nu+p.hopf.ilam)]; 
opt=optimset(p.fsol.opt,'TolFun',p.fsol.tol,'MaxIter',p.fsol.imax);
opt=optimset(opt,'Algorithm','Levenberg-Marquardt'); 
%opt=optimset(opt,'Algorithm','trust-region-reflective');
%opt=optimset(opt,'Algorithm','trust-region-dogleg');
opt=optimset(opt,'Display','iter');
%opt=optimset(opt,'ScaleProblem','Jacobian');
[u,r,exitflag,output,fjac] = fsolve(@fsolextf,ua,opt);

y=u(1:nu*tl); T=u(nu*tl+1); lam=u(nu*tl+2); p.u(p.nu+p.hopf.ilam)=u(end); 
res=norm(r,p.sw.norm); 
iter=output.iterations; 
y=reshape(y,na,tl); % convert to tom-input 
A=fjac; jac=fjac(1:nu*tl,1:nu*tl); 