%% 3 Compos, Yang-Dolnik .. JCP, Oct 2002
% large domain, i.e., k_c=0.7
set(0,'DefaultLineLineWidth', 2);
kv=0:0.7:10; kl=length(kv); ev=zeros(3,kl); % Bru 
Du=0.01; Dv=0.1; Dw=1; c=1; d=1; a=0.95; bv=2.7:0.05:3; bl=length(bv); 
for j=1:bl
b=bv(j);   
for i=1:kl; k=kv(i); 
    A=[[b-1-c-Du*k^2, a^2, d]; 
        [-b, -a^2-Dv*k^2, 0]; 
        [c, 0, -c*a/d-Dw*k^2]]; 
    l=eigs(A); [lso,ix]=sort(real(l),'descend'); ev(:,i)=l(ix);
end
[evs,idx]=sort(abs(real(ev(1,:))));
fprintf('b=%g, real=%g, imag=%g, k=%g\n', bv(j),real(ev(1,idx(1))), imag(ev(1,idx(1))), kv(idx(1)));
figure(6); clf; plot(kv,real(ev(1,:)),'*-k', kv,imag(ev(1,:)),'r'); hold on; 
plot(kv,real(ev(2,:)),'--k', kv,imag(ev(2,:)),'--r'); axis([0 10 -1.1 1.1]);
title(['b=' mat2str(b,3)],'Fontsize',p.plot.fs); pause
end
%% smaller domain, k_c=1.4
set(0,'DefaultLineLineWidth', 2);
kv=0:pi/2.8:10; kl=length(kv); ev=zeros(3,kl); % Bru 
Du=0.01; Dv=0.1; Dw=1; c=1; d=1; a=0.95; bv=2.7:0.05:3; bl=length(bv); 
for j=1:bl
b=bv(j);   
for i=1:kl; k=kv(i); 
    A=[[b-1-c-Du*k^2, a^2, d]; 
        [-b, -a^2-Dv*k^2, 0]; 
        [c, 0, -c*a/d-Dw*k^2]]; 
    l=eigs(A); [lso,ix]=sort(real(l),'descend'); ev(:,i)=l(ix);
end
[evs,idx]=sort(abs(real(ev(1,:))));
fprintf('b=%g, real=%g, imag=%g, k=%g\n', bv(j),real(ev(1,idx(1))), imag(ev(1,idx(1))), kv(idx(1)));
figure(6); clf; plot(kv,real(ev(1,:)),'*-k', kv,imag(ev(1,:)),'r'); hold on; 
plot(kv,real(ev(2,:)),'--k', kv,imag(ev(2,:)),'--r'); axis([0 10 -1.1 1.1]);
title(['b=' mat2str(b,3)],'Fontsize',p.plot.fs); pause
end