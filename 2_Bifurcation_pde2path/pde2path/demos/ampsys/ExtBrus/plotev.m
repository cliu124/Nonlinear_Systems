Dx=0.01; Dy=0.1; Dz=1; p.D=[Dx 0 0;  0 Dy 0;  0 0 Dz]; 
a=1.08; b=3.057; p.par=[a b]; p.uh=[a; b/a; a]; u=sym('u',[3 1]); 
Jf=jacobian(f(u,p),u); Js=subs(Jf,u,p.uh); Js=double(Js); 
k=0:0.1:10; ev=zeros(3,size(k,2)); 
for i=1:size(k,2)
  Lk=Js-p.D*k(i)^2; [~,ew]=eig(Lk); ewd=diag(ew);
   [ews, ix]=sort(real(ewd),'descend');
   ev(1,i)=ewd(ix(1)); ev(2,i)=ewd(ix(2)); ev(3,i)=ewd(ix(3)); 
end
figure(1); clf; pl=ev(1,:); hold on; 
plot(k,real(pl),'k',k,imag(pl),'r'); 
pl=ev(2,:); plot(k,real(pl),'.k',k,imag(pl),'.r'); 
pl=ev(3,:); plot(k,real(pl),'--k',k,imag(pl),'--r'); 
plot(k,zeros(1,size(k,2))); 
axis([0 10 -1.5 1.5]); legend('Re','Im'); box on; 
xlabel('|k|'); 