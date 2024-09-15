function [y,T,lam,res,iter,A,p,jac]=fsol_np(p,y1,T1,lam1,ds)
% fsol_np: interface to call fsolve for Hopf in natural parametrization
global pfs; pfs=p; 
nu=p.nu; tl=p.hopf.tl; yv=reshape(y1(1:nu,:),nu*tl,1); 
ua=[yv; T1; p.u(p.nu+p.hopf.ilam)]; 
opt=optimset(p.fsol.opt,'TolFun',p.fsol.tol,'MaxIter',p.fsol.imax);
opt=optimset(opt,'Algorithm','Levenberg-Marquardt'); 
opt=optimset(opt,'Display','iter');
[u,r,exitflag,output,fjac] = fsolve(@fsolf_np,ua,opt);

y=u(1:nu*tl); T=u(nu*tl+1); p.u(p.nu+p.hopf.ilam)=u(end); 
res=norm(r,p.sw.norm); iter=output.iterations; 
y=reshape(y,nu,tl); % convert to tom-input 
A=fjac; jac=fjac(1:nu*tl,1:nu*tl); lam=lam1;  % fixed lam solver 