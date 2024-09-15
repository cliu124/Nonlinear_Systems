function rplot(varargin)
va=varargin; if ischar(va{1}); p=loadp(va{1},va{2}); tit=[va{1} '/' va{2}]; na=2; 
else p=va{1}; na=1; tit=0; end
R=va{na+1}; nr=va{na+2}; fnr=va{na+3}; 
q=p.u(p.nu+3); po=getpte(p); 
nr=p.np; x=po(1, 1:nr); y=po(2, 1:nr); r=sqrt((x.^2+y.^2)/q); 
%[r1,rm]=min(abs(r-R/sqrt(q))); rm=rm(1) 
[r1,ia]=unique(r); 
figure(fnr); clf; % plot(r(1:rm),p.u(p.np+1:p.np+rm)); 
%plot(r,p.u(1:p.np),'-'); 
plot(r1,p.u(ia),'-'); 
grid on; 
axis tight; set(gca,'fontsize',14); if tit~=0; title(tit); end; 