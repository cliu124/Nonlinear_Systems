function zdia=vegdiagn(p,wnr,pfak,fn)
sol=p.cp; s1=p.oc.s1; u0=p.u0; u1=p.oc.u1; rho=s1.u(s1.nu+1); 
tv=sol.par(1)*sol.t; y=sol.u; sl=length(tv); n=s1.np; zdia=zeros(4,sl); 
for i=1:sl    
    zdia(1,i)=norm(y(1:2*n,i)-u1(1:2*n),'inf'); % diff of states from target 
    zdia(2,i)=jca(s1,y(:,i)); 
    zdia(3,i)=zdia(2,i)*exp(-rho*tv(i)); 
end
plot(tv,zdia(1,:),'k',tv,pfak*zdia(2,:),'m',tv,pfak*zdia(3,:),'g','LineWidth',2);
legend('||v(t)-v_1||_\infty',[mat2str(pfak,2) 'J_{c,a}(t)'],[mat2str(pfak,3) 'e^{-\rho t}J_{c,a}(t)']); 
jp=jcaiT(s1,sol,rho)+disjcaT(s1,sol,rho,0); % jca of path
j0=jca(s1,u0)/rho; j1=jca(s1,s1.u)/rho; 
tit=[fn.sd0 '/' fn.sp0 ' to ' fn.sd1 '/' fn.sp1 '_{.}']; 
fprintf('J0=%g, Jp=%g, J1=%g\n',j0,jp,j1); 
tit={tit, ['J_0=' mat2str(j0,5) ',  J=' mat2str(real(jp),5) ',  J_1=' mat2str(j1,5)]}; 
figure(wnr); set(gca,'FontSize',14); %s1.plot.fs); 
t0=min(tv); t1=max(tv); y0=0; y1=max(max(zdia(1,:)),max(zdia(1,:)))*1.05; 
axis([t0 t1 y0 y1]); xlabel('t');title(tit);
 