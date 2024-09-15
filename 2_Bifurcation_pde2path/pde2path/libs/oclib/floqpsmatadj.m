function [Fu1,Fu2,d]=floqpsmatadj(pd,varargin)
% floqpsmatNT_adj: compute multipliers by pqzSchur; here also return unstable
% eigenspace in Fu2, uses pqzschur-reorder
if(ischar(pd));  if (nargin>1 && ischar(varargin{1})); p=loadp(pd,varargin{1}); else; p=loadp(pd); end
else; p=pd;
end
y=p.hopf.y; tl=p.hopf.tl; nu=p.nu; m=tl-1; t=p.hopf.t; T=p.hopf.T; h=t(2:end)-t(1:end-1);
lam=getlam(p); [~,jac,~,~]=tomassempbc(p,y,T,lam);
AA=zeros(nu,nu,m); BB=AA;
BB(:,:,1)=h(1)*jac(1:nu,1:nu);
AA(:,:,1)=-h(1)*jac(1:nu,(m-1)*nu+1:m*nu);
for i=2:m % rows of jac
    BB(:,:,i)=h(i)*jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag
    AA(:,:,i)=-h(i)*jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag
end
tic; [~,~,~,F1,F2,d] = pqzschur_reorder_adj(AA,BB); toc; % resorted to have unstable F in 1:m/2
Fu2=F2(1:nu/2+1,:); 
Fu1=F1(1:nu/2,:);