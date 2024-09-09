%% relate fold to weierstrass-repr
t=linspace(0,1,200); 
%%
f1=@(t,a) 2*(1-t.^2)./sqrt(t.^8+a*t.^4+1)+4./sqrt(16*t.^4-16*t.^2+2+a); 
f2=@(t,a) 8*t./sqrt(t.^8+a*t.^4+1); 
%%
a=14; y1=f1(t,a); y2=f2(t,a); %mclf(1); plot(t,f1(t,a),t,f2(t,a)); 
E=trapz(t,y1); F=trapz(t,y2); E/F
%%
del1=5.9455; del2=6.399; al=5; ar=50;
a=linspace(al,ar,100); del=0*a; 
for i=1:length(a); 
    y1=f1(t,a(i)); y2=f2(t,a(i));
    E(i)=trapz(t,y1); F(i)=trapz(t,y2); del(i)=2*pi*F(i)/E(i); 
end
%%
mclf(2); plot(a,del,'b'); hold on; plot(a, del1+0*a,'m'); 
plot(a,del2+0*a,'r'); a1=7.402; a2=28.78; fs=14; 
plot([a1, a1],[5 7],'k'); plot([a2, a2],[5 7],'k'); 
axis([al ar min(del) 1.01*max(del)]); 
grid on; legend('2\pi F/E','\delta_1','\delta_f'); 
xlabel('a'); ylabel('\delta'); xticks([20 40]); 
text(a1,0.99*min(del),'a_1','fontsize',16); 
text(a2,0.99*min(del),'a_2','fontsize',16); 
set(gca,'fontsize',16);