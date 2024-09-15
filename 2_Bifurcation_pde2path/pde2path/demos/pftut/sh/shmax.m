muc=0;muf=-0.2;
mumu=linspace(muf,muc,50);
EEp=zeros(1,50);
a=3;
b=6;
for i=1:50
mu=mumu(i);
As=1/15-sqrt(1/15^2+mu/7.5/6)
Ep=mu/2*6*As^2+0.5*a*4*As^3+b/8*6*As^4+b/2*15*As^4+3*b*As^4;
EEp(i)=Ep;
end
%%
figure(1);clf;
plot(mumu,real(EEp),'linewidth',2);hold on;
plot([mumu(1),mumu(end)],[0,0],'color','k');
axis([real(mumu(1)),real(mumu(end)),real(EEp(1)),real(EEp(end))]);
set(gca,'fontsize',25)