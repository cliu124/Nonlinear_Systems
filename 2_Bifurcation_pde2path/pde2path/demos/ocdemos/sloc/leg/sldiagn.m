function dia=sldiagn(sol,wnr) % standard diagnostics, plot u_1, lam_1
% use as template to plot adapted diagnostics 
global s0 u0 u1; 
x=sol.x; y=sol.y; sl=length(x); np=s0.np; dia=zeros(2,sl); 
for i=1:sl
    dia(1,i)=norm(y(1:np,i)-u0(1:np),'inf'); % diff state from IC 
    dia(2,i)=norm(y(1:np,i)-u1(1:np),'inf'); % diff state from target
    dia(3,i)=norm(y(np+1:2*np,i)-u0(np+1:2*np),'inf'); % diff costate from IC
    dia(4,i)=norm(y(np+1:2*np,i)-u1(np+1:2*np),'inf'); % diff costate from target
end
figure(wnr); set(gca,'FontSize',s0.plot.fs); 
plot(x,dia(1,:),'k',x,dia(2,:),'b',x,dia(3,:),'m',x,dia(4,:),'g','LineWidth',2); 
legend('||P-P_0||_\infty','||P-P_1||_\infty','||q-q_0||_\infty','||q-q_1||_\infty'); 
axis tight;
end 