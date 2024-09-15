function w=asT(x,y,z,p,ss)  
% asT: trilinear form used in ampsys 
u=sym('u',[ss 1]);w=0;
if isfield(p,'sb'); sb=p.sb; else; sb=0; end
for i=1:ss
    for j=1:ss
        for k=1:ss
            if i<=j && j<=k
                df=diff(f(u,p),u(i));df=diff(df,u(j));df=diff(df,u(k));df=subs(df,u,p.uh);
                if sb==0;df=double(df);end
            end
            if i==j && j==k
                w=w+1/6*df*(x(i)*y(j)*z(k));
            end
            if i<j && j<k
               o=x(i)*y(j)*z(k)+x(i)*z(j)*y(k)+...
                   y(i)*x(j)*z(k)+y(i)*z(j)*x(k)+...
                   z(i)*x(j)*y(k)+z(i)*y(j)*x(k);
               w=w+1/6*df*o;
            end
            if (i==j && j~=k) || (i~=j && j==k)
                o=x(i)*y(j)*z(k)+x(i)*z(j)*y(k)+...
                   y(i)*x(j)*z(k)+y(i)*z(j)*x(k)+...
                   z(i)*x(j)*y(k)+z(i)*y(j)*x(k);
               w=w+0.5*1/6*df*o;
            end                          
        end
    end
end