function  mov=homovplot(p,cnr,wnr,st) 
z=p.hopf.y; tv=p.hopf.t; T=p.hopf.T; tl=length(tv); %figure(wnr); clf
zp=z(:,:); movc=1; 
for i=1:tl
    plotsolu(p,zp(:,i),wnr,cnr,st); view(20,60); colormap cool; 
    set(gca,'XTick',[]); set(gca,'yTick',[]); set(gca,'zTick',[]);
    zlabel(''); fprintf('t=%g,  ', T*tv(i)); title(['t=' mat2str( T*tv(i),3)]); 
    mov(movc)=getframe(wnr); movc=movc+1; 
    pause
end
    