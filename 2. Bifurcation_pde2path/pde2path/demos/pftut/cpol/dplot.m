k=0.05; ga=1; u=0:0.05:2.4; w=1.6:0.05:2.6; [U,W]=meshgrid(u,w);
f=(k+ga*U.^2./(1+U.^2)).*W-U;
%
figure(1); clf; 
[M,c]=contour(U,W,f,[0 0],'ShowText','off'); %,[0 0]); 
c.LineWidth=2; 
xlabel('u'); ylabel('w'); legend('f(u,w)=0'); 
%%
u=0:0.05:1.8; w=2; 
f=(k+ga*u.^2./(1+u.^2)).*w-u;
figure(2); clf; plot(u,f); xlabel('u'); ylabel('f(u,2)'); axis tight;