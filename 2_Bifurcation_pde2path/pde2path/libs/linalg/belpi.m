function [x,y,p,r]=belpi(A,b,c,d,f,g,p)
% BELPI: bordered elimination with post iterations
%
%  [x,y,p]=belpi(A,b,c,d,f,g,p)
%
% called by lssbel, calls bel 
% see also: lss, blss, lssbel, blssbel, bel
[x,y,p]=bel(A,b,c,d,f,g,p); % block-elim  
fs=f-(A*x+b*y); gs=g-(c*x+d*y);
r=norm([fs;gs],'inf'); iter=0;rs=r;
while r>p.bel.tol && iter<p.bel.imax
[x1,y1,p]=bel(A,b,c,d,fs,gs,p); x=x+x1;y=y+y1; iter=iter+1; % corrector
fs=f-(A*x+b*y); gs=g-(c*x+d*y); r=norm([fs;gs],'inf');
rs=[rs,r]; 
end
if iter>0 && p.sw.verb>1
   str1=[num2str(iter),' post iterations done in belpi'];
   str2=['residuals were ',num2str(rs)];
   disp(str1);disp(str2);
end

