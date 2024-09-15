function p=doublesec(p,q,phi0) % ad_hoc function to double a sector: take sol 
% from q, rotate by -phi0, then reflect at y=0 and put into p
np1=q.np; uq=q.u(1:np1); % # of points in q, and solution  
po1=getpte(q); x1=po1(1,:); y1=po1(2,:); % points in q
x1a=cos(phi0)*x1+sin(phi0)*y1; y1a=-sin(phi0)*x1+cos(phi0)*y1; % rotate by -phi0
po2=getpte(p); x2=po2(1,:); y2=po2(2,:); np2=p.np; % points in p
if phi0>0; i1=find(y2<0); else i1=find(y2>0); end % find points of p in upper/lower half 
i2=setdiff(1:np,i1); x2a=x2(i1); y2a=y2(i1); % half the points in p

u2a=p2interpol(x2(i1),y2(i1),uq,x1a,y1a,q); % interpolate q to p (one half) 
 
i2l=length(i2); u2b=zeros(i2l,1); % fill the other half by reflection at y=0 
for j=1:i2l; % probably not very efficient ... 
    i=i2(j); x=x2(i); y=y2(i); % point cordinates 
    [m,in]=min((x2a-x).^2+(y2a+y).^2); % index of mirror image 
    u2b(j)=-u2a(in);  % reflect soln 
end
u2=zeros(np2,1); u2(i1)=u2a; u2(i2)=u2b; % fill soln in p 
p.u(1:p.nu)=u2; p.u(p.nu+1:end)=q.u(q.nu+1:end); % take pars from q and check resi 
r=pderesi(p,p.u); p=resetc(p); p.sw.verb=2; fprintf('initial residual=%g\n', norm(r,'inf')); 
% correct by a Newton loop 
[u,res,iter,Gu,Glam,p]=nloop(p,p.u);  p.u(1:p.nu)=u(1:p.nu); 
fprintf('residual after Newton=%g\n',res); 


