%%
% floqpqzsap: compute multipliers by pqzSchur
%  here adapted to pollution-OC: log-plots, and d=nu/2-ind; 
%
%  [muv1, muv2,ind,h]=floqpqzsap(pd,varargin)
function [muv1, muv2, ind,h]=floqpsap(pd,varargin)
if(ischar(pd)) 
    if (nargin>1 && ischar(varargin{1})); p=loadp(pd,varargin{1});
    else p=loadp(pd); end
else p=pd; end
y=p.hopf.y; tl=p.hopf.tl; nu=p.nu; m=tl-1; %na=m*tl; M=p.mat.M; 
t=p.hopf.t; T=p.hopf.T; h=t(2:end)-t(1:end-1); %par=p.u(p.nu+1:end); 
lam=getlam(p);
[f,jac,f_T,f_lam]=tomassempbc(p,y,T,lam); %figure(1); spy(jac)
AA=zeros(nu,nu,m); BB=AA; 
BB(:,:,1)=h(1)*jac(1:nu,1:nu);
AA(:,:,1)=-h(1)*jac(1:nu,(m-1)*nu+1:m*nu); 
for i=2:m; % rows of jac
  BB(:,:,i)=h(i)*jac((i-1)*nu+1:i*nu, (i-1)*nu+1:i*nu); % diag 
  AA(:,:,i)=-h(i)*jac((i-1)*nu+1:i*nu, (i-2)*nu+1:(i-1)*nu); % subdiag
end 
tic; [A,B,Q,Z] = pqzschur(AA,BB); toc
d=ones(nu,1); 
for i=1:nu
  for k=1:m; 
    d(i)=d(i)*A(i,i,k)/B(i,i,k);
  end
end
nprint=4; 
[ds,idx]=sort(abs(d)); d=d(idx); 
for i=1:nprint; fprintf('%g+%gi,    ',real(d(i)),imag(d(i))); end; fprintf('\n'); 
[mu1,idtr]=min(abs(abs(d)-1)); mu1=d(idtr); 
muv1=d(1:idtr); muv2=d(idtr+1:end);  
try
for i=1:nprint; fprintf('%g+%gi,    ',real(muv1(end-i)),imag(muv1(end-i))); end; fprintf('\n'); 
for i=1:nprint; fprintf('%g+%gi,    ',real(muv2(i)),imag(muv2(i))); end; fprintf('\n'); 
end
figure(7); clf; set(gca,'FontSize',12); 
ind=nu/2-length(muv1); 
h=plot(real(muv1), imag(muv1),'*b'); 
set(h, 'Markersize',8); hold on;
if ind>0; h=plot(real(muv2(1:ind)),imag(muv2(1:ind)),'*r'); set(h, 'Markersize',8); end
s=0:0.05:2*pi; plot(cos(s),sin(s)); axis image;
figure(8); clf;  set(gca,'FontSize',12); 
h=semilogy(1:length(muv2),abs(muv2),'-*'); axis tight; set(h, 'Markersize',8);
fprintf('error in ga_1=%g\n', mu1-1); 
end
