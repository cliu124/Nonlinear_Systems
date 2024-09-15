% for OC-Hopf 
% rem: bradat=[count; ptype; ineg; getlam; err; |u_1|_L^2] (6 values) 
% here: par(1..5), T, min, max, |u(.,.)|_L^2, J_c(CSS)/J_c(u_H), J_c(u_H(.+T/2))
%                  0                                                  0             
function out=pollbra(p,u) 
gotit=0; nu=p.nu; par=u(p.nu+1:end); rho=par(1); 
if isfield(p,'hopf'); ho=p.hopf; 
  if isfield(ho,'y'); T=ho.T; y=ho.y; 
    m1=T; m2=max(max(y(1:nu,:))); m3=min(min(y(1:nu,:))); 
    l2v=zeros(1,ho.tl); jcv=zeros(1,ho.tl); % for num eval of J
    for i=1:length(ho.t) 
        l2v(i)=y(1:nu,i)'*(ho.tom.M(1:nu,1:nu)*y(1:nu,i)); 
        ua=[y(1:nu,i); par]; jcv(i)=jca(p,ua);
    end 
    l2=trapz(ho.t,l2v); m4=sqrt(l2/(p.vol)); 
    tv=T*ho.t(1:length(ho.t)-1); tt=[tv, tv+T, tv+2*T, tv+3*T];   % 4 periods of t
    jcvo=jcv(1:length(ho.t)-1); jcvv=[jcvo,  jcvo, jcvo, jcvo];   % 4 periods of sol 
    s=round(ho.tl/2); jcvs=[jcv(s+1:length(ho.t)-1),jcv(1:s)] ; 
    jcvvs=[jcvs, jcvs, jcvs, jcvs]; 
    set(0,'defaultlinelinewidth',2);
    figure(5); plot(tt,jcvv.*exp(-0*rho*tt),'k'); 
    fs=p.plot.fs; axis([0 T min(jcvv) max(jcvv)]); xlabel('t','fontsize',fs); 
    ylabel('J_{c,a}(t)','fontsize',fs); set(gca,'FontSize',p.plot.fs); 
    set(0,'defaultlinelinewidth',1);
    m5=trapz(tt,exp(-rho*tt).*jcvv); 
    m6=trapz(tt,exp(-rho*tt).*jcvvs); 
    gotit=1; 
  end; 
end 
if ~gotit; m1=0; m2=max(u(1:nu)); m3=min(u(1:nu));
    l2=u(1:nu)'*(p.mat.M*u(1:nu)); m4=sqrt(l2/p.vol); m5=jca(p,u)/rho; m6= jca(p,u)/rho; 
end
out=[u(p.nu+1:p.nu+5); % parameters 
      m1; m2; m3; m4; m5; m6]; 