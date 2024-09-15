function zdia=vegdiagn(sol,wnr,pfak,fn)
global s0 s1 u0 u1; 
x=sol.x; y=sol.y; sl=length(x); n=s1.np; zdia=zeros(4,sl); 
for i=1:sl
    zdia(1,i)=norm(y(1:2*n,i)-u0(1:2*n),'inf'); % diff of states from IC 
    zdia(2,i)=norm(y(1:2*n,i)-u1(1:2*n),'inf'); % diff of states from goal
    zdia(3,i)=jca(s1,y(:,i)); 
    zdia(4,i)=1; % diff of costates from goal 
end
if 1
plot(x,zdia(1,:),'m',x,zdia(2,:),'k',x,pfak*zdia(3,:),'b','LineWidth',2);
legend('||v(t)-v_0||_\infty','||v(t)-v_1||_\infty',[mat2str(pfak,2) 'J_{c,a}']); 
else 
plot(x,zdia(2,:),'k',x,pfak*zdia(3,:),'b','LineWidth',2)
legend('||u(t)-u_1||_\infty',[mat2str(pfak,2) 'J_{c,a}']); 
end 
rho=s1.u(s1.nu+1); jp=jcai(s1,sol,rho)+disjca(s1,sol,rho); % jca of path
j0=jca(s0,s0.u)/rho; j1=jca(s1,s1.u)/rho; 
tit=[fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 '_{.}']; 
fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); 
tit={tit, ['J_0=' mat2str(j0,5) ',  J=' mat2str(real(jp),5) ',  J_1=' mat2str(j1,5)]}; 
figure(wnr); set(gca,'FontSize',14); %s1.plot.fs); 
t0=min(x); t1=max(x); y0=0; y1=max(max(zdia(1,:)),max(zdia(1,:)))*1.05; 
axis([t0 t1 y0 y1]); xlabel('t');title(tit);
 