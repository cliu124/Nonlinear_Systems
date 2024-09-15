function mov=sol2mov(p,sol,wnr)
sl=length(sol.x); movc=1; [po,tr,e]=getpte(p); 
x1=min(po(1,:)); x2=max(po(1,:)); 
y1=min(po(2,:)); y2=max(po(2,:)); 
u1=sol.y(1:p.np,:); u2=-1./sol.y(p.np+1:2*p.np,:); 
z11=min(min(u1)); z12=max(max(u1)); z21=min(min(u2)); z22=max(max(u2));
ax1=[x1 x2 y1 y2 z11 z12]; ax2=[x1 x2 y1 y2 z21 z22]; 
for i=1:sl; 
    plotsolusf(p,sol.y(:,i),wnr,1,p.plot.pstyle,13,1,[2 3 4 5]); % vert
    title(['t=' mat2str(sol.x(i),3)]); 
    xlabel(''); ylabel(''); zlabel(''); 
    text(x1-4,y1,z12,'P','FontSize',16); 
    axis(ax1); view(14,60); freezeColors; 
    u=-1./sol.y(:,i); 
    plotsolusf(p,u,wnr,2,p.plot.pstyle,13,1,[9 10 11 12]); 
    xlabel(''); ylabel(''); zlabel('v');  colormap hot; 
    text(x1-4,y1,1.2*z22,'k','FontSize',16); 
    axis(ax2); view(14,60); pause
    mov(movc)=getframe(wnr);movc=movc+1; 
end
%movie2avi(mov,'m1.avi');
% on linux: convert with 
% mencoder in.avi -o out.avi -speed 0.2 -ovc lavc 
% or 
% ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% to about 1/100 in sizega=p.u(p.nu+4); u=p.u; n=p.np; 