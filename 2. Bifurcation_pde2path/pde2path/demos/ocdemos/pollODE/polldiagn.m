function dia=polldiagn(p,wnr,varargin) 
% use as template to plot adapted diagnostics 
f1=1; f2=1; f3=1; f4=1; 
try; f1=varargin{1}; f2=varargin{2}; f3=varargin{3}; f4=varargin{4}; catch; end; 
s1=p.oc.s1; u0=p.oc.u0; u1=p.oc.u1;
sol=p.cp; T=sol.par(1); 
t=T*sol.t; y=sol.u; sl=length(t); n=s1.nu/4; dia=zeros(8,sl); 
for i=1:sl
    dia(1,i)=norm(y(1:n,i)-u0(1:n),'inf'); % diff v from init 
    dia(2,i)=norm(y(1:n,i)-u1(1:n),'inf'); % diff v from target
    dia(3,i)=norm(y(n+1:2*n,i)-u0(n+1:2*n),'inf'); % diff w from IC
    dia(4,i)=norm(y(n+1:2*n,i)-u1(n+1:2*n),'inf'); % diff w from target  
    dia(5,i)=norm(y(2*n+1:3*n,i)-u1(2*n+1:3*n),'inf'); % diff lam from target  
    dia(6,i)=norm(y(3*n+1:4*n,i)-u1(3*n+1:4*n),'inf'); % diff mu from target   
end
figure(wnr); set(gca,'FontSize',s1.plot.fs); 
plot(t,f1*dia(2,:),'k',t,f2*dia(4,:),'b',t,f3*dia(5,:),'m',t,f4*dia(6,:),'g','LineWidth',2); 
sw=2; 
switch sw
    case 1; legend([mat2str(f1) '||v-v_1||_\infty'],[mat2str(f2) '||w-w_1||_\infty'],...
    [mat2str(f3) '*||\lambda_1-\lambda_1||_\infty'],[mat2str(f4) '*||\mu-\mu_2||_\infty']);   
    case 2; legend([mat2str(f1) '||v-v_1||_\infty'],[mat2str(f2) '||w-w_1||_\infty'],...
    ['||\lambda_1-\lambda_1||_\infty'],['||\mu-\mu_2||_\infty']);   
end
axis tight; xlabel('t'); set(gca,'fontsize',14); 
rho=s1.u(s1.nu+p.oc.rhoi); 
[Jtran,jcv,jcvd]=jcaiT(s1,sol,rho); 
figure(12); clf; 
 plot(t,jcv,'k'); axis tight; legend('J_{ca}'); %,'e^{-\rho t}J_{ca}'); 
 title(['J=' mat2str(Jtran,4)]); 
axis tight; xlabel('t'); set(gca,'fontsize',14); 
end 