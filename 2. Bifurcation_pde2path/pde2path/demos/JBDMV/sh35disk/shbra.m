function out=shbra(p,u)
M=getM(p); n=p.np; J=shJ(p,u); po=getpte(p);
v=p.u(1:p.nu/2); radr=find(abs(po(2,:))<1e-2); xy=[po(1,radr); v(radr)'];
xy=xy'; xy=sortrows(xy,1); xmax=max(xy(:,1)); 
%figure(10); clf; plot(xy(:,1),xy(:,2),'-')
out=[u(p.nu+1:end);
    sqrt(trapz(xy(:,1),xy(:,2).^2))/sqrt(xmax); 
     max(v);
     min(v);
     J; sqrt(u(1:n)'*(M(1:n,1:n)*u(1:n)))/sqrt(p.Om)];
     %sqrt(trapz(xy(:,1),xy(:,2).^2/30))]; 
    
% out=[max(v);sqrt(trapz(xy(:,1),xy(:,2).^2/30))];