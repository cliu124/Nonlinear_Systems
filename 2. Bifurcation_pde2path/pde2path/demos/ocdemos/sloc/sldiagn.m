function dia=sldiagn(p,wnr) % diagnostics for SLOC model, new syntax 
% use as template to plot adapted diagnostics 
s1=p.oc.s1; u1=p.oc.u1; sol=p.cp; T=sol.par(1); 
t=T*sol.t; y=sol.u; tl=length(t); np=s1.np; dia=zeros(2,tl); 
for i=1:tl
    dia(1,i)=norm(y(1:np,i)-u1(1:np),'inf'); % diff state from IC 
    dia(2,i)=norm(y(np+1:2*np,i)-u1(np+1:2*np),'inf'); % diff costate from target   
end
rho=s1.u(s1.nu+1); [jcaval,jcav,jcavd]=jcaiT(s1,sol,rho); 
jp=jcaval+disjcaT(s1,sol,rho,0); % jca of path
j0=jca(s1,p.u0)/rho; j1=jca(s1,s1.u)/rho; 
%tit=[fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 '_{.}']; 
figure(wnr); set(gca,'FontSize',s1.plot.fs); 
psw=1; 
switch psw; 
  case 1; plot(t,dia(1,:),'k',t,jcav(:),'m',t,jcavd(:),'g','LineWidth',2); 
  legend('||v-v_1||_\infty','J_{c,a}','e^{-\rho t}J_{c,a}'); 
  case 3; plot(t,dia(1,:),'k',t,dia(2,:),'b',t,jcav(:),'m',t,jcavd(:),'g','LineWidth',2); 
  legend('||v-v_1||_\infty','||\lambda-\lambda_1||_\infty','J_{c,a}','e^{-\rho t}J_{c,a}'); 
end 
xlabel('t'); fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); 
tit={['J_0=' mat2str(j0,3) ',  J=' mat2str(real(jp),3) ',  J_1=' mat2str(j1,3)]}; 
figure(wnr); set(gca,'FontSize',14); title(tit); axis tight;
