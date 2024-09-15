function y0d=sety0dot(p,y,par,T)
% sety0dot: set y0dot for the phase condiion (in t) for Hopf orbits
n=p.np; tl=p.hopf.tl; M=p.mat.M(1:2*n,1:2*n); 
yd1=zeros(p.nu,p.hopf.tl); yd2=yd1; y0d=yd1; 
y0dsw=2; if isfield(p.hopf,'y0dsw'); y0dsw=p.hopf.y0dsw; end 
scG=1; if isfield(p.hopf,'scG'); scG=p.hopf.scG; end 
switch y0dsw
    case 0; % use PDE (consistent, if (at least first compo) is evol. eqn)        
      for i=1:tl;  p.i=i; f=hosrhs(0,y(:,i),p,par,T); %y0d(:,i)=scG*f(:); 
      nn=n; y0d(1:nn,i)=scG*M(1:nn,1:nn)\f(1:nn);
      nn=n; y0d(2*n+1:3*n,i)=scG*M(1:nn,1:nn)\f(2*n+1:3*n);
      end;     
    case 1; % 1st order FD at t=0,1        
      i1=1; i2=3*n+1; iv=[i1:i2]; 
      h=p.hopf.t(2:end)-p.hopf.t(1:end-1);   % just forward FD 
      for i=1:tl-1; yd1(iv,i)=(y(iv,i+1)-y(iv,i))./h(i); end
      for i=2:tl; yd2(iv,i)=(y(iv,i)-y(iv,i-1))./h(i-1); end
      y0d=0.5*(yd1+yd2); % 2nd-order in the middle
      y0d(:,1)=yd1(:,1); 
      y0d(:,end)=yd2(:,end); 
    case 2; % 2nd order FD         
      h=p.hopf.t(2:end)-p.hopf.t(1:end-1); h=[h h(end)]; 
      cmp=1; n1=(cmp-1)*n+1; n2=cmp*n; 
      for i=2:tl-1; y0d(n1:n2,i)=(y(n1:n2,i+1)-y(n1:n2,i-1))./(2*h(i)); end
      %  y0d(1:n,1)=(y(1:n,2)-y(1:n,end))./h(1);
      y0d(n1:n2,1)=(y(n1:n2,2)-y(n1:n2,1))./h(1); % just forward
      y0d(n1:n2,end)=(y(n1:n2,end)-y(n1:n2,end-1))./h(end); 
       %y0d(:,end)=y0d(:,1); 
     case 3; % via makeDt 
        t=p.hopf.t; 
        Dt=makeDtm(t); i1=1; i2=n; 
        for i=i1:i2; y0d(i,:)=(Dt*y(i,:)')'; end;     
end
