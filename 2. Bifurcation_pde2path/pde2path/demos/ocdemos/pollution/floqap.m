%%
% floqap: a posteriori compute multipliers and instability index of u_H 
%  
% compute Floquet multipliers and instability index of u_H 
% in p (p2p-struct, or directory) 
% muv1,muv2=computed multipliers, ind=ind(u_H)
% mo=monodromy, h=plot-handle (to modify marker-size etc) 
function [muv1, muv2,ind,mo,h]=floqap(pd,varargin)
t1=tic; noa=nargin-1; k=1; 
if(ischar(pd)) 
    if (noa>0 && ischar(varargin{1})); p=loadp(pd,varargin{1}); k=k+1; noa=noa-1; 
    else p=loadp(pd); 
    end
else p=pd; end
n1=p.nu; if noa>0; n1=varargin{k}; end  % Nr of Multipl. to compute
lam=getlam(p); 
y=p.hopf.y; tl=p.hopf.tl; nu=p.nu; m=tl-1; 
t=p.hopf.t; T=p.hopf.T; h=t(2:end)-t(1:end-1); 
[f,jac]=tomassempbc(p,y,T,lam);
AA=zeros(nu,nu,m); BB=AA; 
BB(:,:,1)=h(1)*jac(1:nu,1:nu);
AA(:,:,1)=-h(1)*jac(1:nu,(m-1)*nu+1:m*nu); 
for i=2:m; % rows of jac
  BB(:,:,i)=h(i)*jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag 
  AA(:,:,i)=-h(i)*jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag
end 
mo=speye(nu); for k=1:m; fak=BB(:,:,k)\AA(:,:,k); mo=fak*mo; end; % building M 
ok=0; 
while ~ok; 
  [V,mu,flag]=eigs(mo,n1,'lm');  %  multipliers 
  muv=mu*ones(1,n1)'; [muvs,ix]=sort(abs(muv)); muv=muv(ix); ok=1; 
  nprint=min(6,n1); 
  for i=1:nprint; fprintf('%g+%gi,    ',real(muv(i)),imag(muv(i))); end; fprintf('\n');
end
muv(1)
idx=find(abs(muv)>1+p.hopf.fltol); 
muv2=muv(idx); muv1=setdiff(muv,muv2); 
lm=length(muv1); mu1=muv1(end); 
figure(7); clf; set(gca,'FontSize',12); 
h=plot(real(muv1), imag(muv1),'*b'); set(h, 'Markersize',8); hold on;
s=0:0.1:2*pi; plot(cos(s),sin(s)); axis image;
figure(8); clf;  set(gca,'FontSize',12); 
h=semilogy(1:length(muv2),abs(muv2),'-*'); axis tight; set(h, 'Markersize',8);
ind=length(muv2); 
t1=toc(t1); 
fprintf('floqap: %g seconds, error in mu_1=%g, #unst=%i\n', t1, mu1-1, ind); 
end
 