function [muv1,muv2,ind,h]=floqpsap(pd,varargin)
% floqpsap: Floquet multipliers by pqzSchur (a posteriori) 
%
%  [muv1, muv2, ind, h]=floqpsap(pd,varargin)
% pd=(p2p-struct, or directory) 
% muv1,muv2=computed multipliers, ind=ind(u_H)
% h=plot-handle (to modify marker-size etc) 
%
noa=nargin-1; k=1; 
if(ischar(pd)) 
    if (nargin>1 && ischar(varargin{1})); p=loadp(pd,varargin{1}); 
        k=k+1; noa=noa-1;
    else p=loadp(pd); end
else p=pd; 
end
aux=[]; if noa>0; aux=varargin{k}; end 
wnr=7; if isfield(aux,'wnr'); wnr=aux.wnr; end 
t1=tic; tl=p.hopf.tl; nu=p.nu; m=tl-1; 
lam=getlam(p); [f,jac,f_T,f_lam]=tomassempbc(p,p.hopf.y,p.hopf.T,lam); 
[AA,BB]=flmcoll(p,jac); 
tic; [A,B,Q,Z] = pqzschur(AA,BB); toc
d=ones(nu,1); 
for i=1:nu
  for k=1:m; 
    d(i)=d(i)*A(i,i,k)/B(i,i,k);
  end
end
[ds,idx]=sort(abs(d)); d=d(idx); 
idx=find(abs(d)>1+p.hopf.fltol); muv2=d(idx); % unstable multipliers
muv1=setdiff(d,muv2); % the stable multipliers 
ind=length(muv2); 
[m1,idx]=min(abs(muv1-1)); mu1=muv1(idx); mu1err=mu1-1;  % should be 1
t1=toc(t1); 
fprintf('floqps, %g seconds, error in ga_1=%g, #unst=%i\n', t1, mu1err, ind); 
if abs(mu1err)>p.hopf.fltol; 
    fprintf('floqps: error in ga1 is probably too large ..\n'); end
figure(wnr); clf; set(gca,'FontSize',12); 
h=plot(real(muv1), imag(muv1),'*b', real(muv2), imag(muv2), '*r'); 
set(h, 'Markersize',8);
hold on; s=linspace(0,2*pi,100); plot(cos(s),sin(s)); axis image; %pause
end
