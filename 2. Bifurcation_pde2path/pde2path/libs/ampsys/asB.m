function w=asB(x,y,p,ss) 
% asB: bilinear form used in ampsys
u=sym('u',[ss 1]);w=0;
if isfield(p,'sb'); sb=p.sb; else; sb=0; end
for i=1:ss
    for j=1:ss
        if i<=j
        df=diff(f(u,p),u(i));df=diff(df,u(j)); % \pa_{ij} f
        df=subs(df,u,p.uh);
        if sb==0;df=double(df);end
        end
    if i==j
    w=w+0.5*df*x(i)*y(j);
    else
    if i<j
       w=w+0.5*df*(x(i)*y(j)+x(j)*y(i)); 
    end 
    end
    end
end
