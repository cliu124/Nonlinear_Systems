function zdia=lvdiagn(sol,wnr,fn,sw) % diagnostics for lv-bd-control 
% sw=1: J,  2: k_1,2,  3: v_1,2,   4: lam_1,lam_2,  5: h1,h2
global s0 s1; 
set(0,'defaultlinelinewidth',2)
x=sol.x; y=sol.y; sl=length(x); n=s1.np; 
switch sw; case 1; zdia=zeros(1,sl); otherwise; zdia=zeros(2,sl); end 
par=s1.u(s1.nu+1:end); al1=par(9); al2=par(10); 
p1=par(11); p2=par(12); c1=par(13); c2=par(14); 
for i=1:sl
  switch sw
      case 1; zdia(1,i)=jca(s1,[y(:,i); par]); 
      case 2; k=cficon(s1,[y(:,i); par]); zdia(1,i)=k(1); zdia(2,i)=k(2); 
      case 3; zdia(1,i)=y(1,i); zdia(2,i)=y(n+1,i); 
      case 4; zdia(1,i)=y(2*n+1,i); zdia(2,i)=y(3*n+1,i);    
      case 5; k=cficon(s1,[y(:,i); par]); v1=y(1,i); v2=y(n+1,i); 
          zdia(1,i)=v1^al1*k(1)^(1-al1); 
          zdia(2,i)=v2^al2*k(2)^(1-al2); 
  end
end
figure(wnr); 
switch sw
  case 1; loglog(x,zdia(1,:),'b'); legend('J_c'); 
      %semilogy(x,zdia(1,:),'b'); legend('J_c'); 
      %plot(x,zdia(1,:),'b'); legend('J_c'); 
  case 2; %plot(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('k_1','k_2'); 
      loglog(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('k_1','k_2'); 
  case 3; %plot(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('v_1','v_2'); 
      loglog(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('v_1','v_2'); 
  case 4; plot(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('\lambda_1','\lambda_2');  
      %loglog(x,zdia(1,:),'b',x,zdia(2,:),'r'); legend('\lambda_1','\lambda_2');  
  case 5; plot(x,zdia(1,:),'b',x,zdia(2,:),'r',x,0*zdia(1,:),'k--'); legend('h_1','h_2');  
      %loglog(x,zdia(1,:),'b',x,zdia(2,:),'r',x,0*zdia(1,:),'k--'); legend('h_1','h_2'); 
end
rho=s1.u(s1.nu+8); jp=jcai(s1,sol,rho)+disjca(s1,sol,rho); % jca of path
j0=jca(s0,s0.u)/rho; j1=jca(s1,s1.u)/rho; 
tit=['Eq to ' fn.sd1 '/' fn.sp1 ', J=' mat2str(real(jp),4)]; 
fprintf('Jp=%g, J1=%g\n',jp,j1); 
t0=min(x); t1=max(x); 
y0=min(min(zdia(:,:))); y1=max(max(zdia(:,:)))*1.05; 
if y0>0; y0=0.95*y0; else y0=1.05*y0; end; axis([t0 t1+0.1 y0 y1]); 
xlabel('t'); %title(tit); 
set(gca,'FontSize',s1.plot.fs); set(0,'defaultlinelinewidth',1)