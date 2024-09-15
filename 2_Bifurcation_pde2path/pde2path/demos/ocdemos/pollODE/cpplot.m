function cpplot(poc,wnr,varargin)
f1=1; f2=10; 
try; f1=varargin{1}; f2=varargin{2}; f3=varargin{3}; f4=varargin{4}; catch; end; 
ts=1; if nargin>2; ts=varargin{1}; end 
cp=poc.cp; v=[70,30]; s1=poc.oc.s1; 
np=s1.np; T=real(cp.par); tv=T*cp.t(1:ts:end); tl=length(tv); 
z=cp.u; par=s1.u(s1.nu+1:end); ga=par(5); rho=par(1); 
jcv=zeros(1,tl); kv=jcv; l1v=jcv; l2v=jcv; 
for i=1:tl;
    u=[z(:,i); par]; jc=polljcf(s1, u); jcv(i)=jc(1); 
    l1=z(3,i); kv(i)=-(1+l1)/ga;
    l2=z(4,i); l1v(i)=l1; l2v(i)=l2; 
end
figure(wnr); clf 
plot(tv, jcv,'-r', tv,z(1,:),'-k',tv,z(2,:),'-b',tv,20*kv,'-m'); 
legend('J_c','v_1','v_2','20*q');  xlabel('t'); 
axis tight; [Jtran,jcv,jcvd]=jcaiT(s1,cp,rho); 
title(['J=' mat2str(Jtran,4)]); set(gca,'fontsize',14); 
figure(11); clf 
plot(tv, jcv,'-r', tv,z(1,:),'-k',tv,z(2,:),'-b',tv,20*kv,'-m'); xlabel('t'); 
legend('J_c','v_1','v_2','20*q'); % axis([0 75 -0.3 1]); set(gca,'fontsize',14); 
