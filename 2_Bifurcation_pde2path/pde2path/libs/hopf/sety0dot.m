function y0d=sety0dot(p,y,par,T)
% sety0dot: set y0dot for the phase condiion (in t) for Hopf orbits
n=p.nu/p.nc.neq; tl=p.hopf.tl; 
yd1=zeros(p.nu,p.hopf.tl); yd2=yd1; y0d=yd1; 
y0dsw=2; if isfield(p.hopf,'y0dsw'); y0dsw=p.hopf.y0dsw; end 
scG=1; if isfield(p.hopf,'scG'); scG=p.hopf.scG; end 
switch y0dsw
    case 0; % use PDE (consistent, if (at least first compo) is evol. eqn) 
      for i=1:tl-1;  f=hosrhs(0,y(:,i),p,par,T); %y0d(:,i)=scG*f(:); 
      y0d(1:n,i)=scG*f(1:n);
      end; 
    case 1; % 1st order FD at t=0,1
      h=p.hopf.t(2:end)-p.hopf.t(1:end-1);   % just forward FD 
      for i=1:tl-1; yd1(1:n,i)=(y(1:n,i+1)-y(1:n,i))./h(i); end
      for i=2:tl; yd2(1:n,i)=(y(1:n,i)-y(1:n,i-1))./h(i-1); end
      y0d=0.5*(yd1+yd2); % 2nd-order in the middle
    case 2; % 2nd order FD 
      h=p.hopf.t(2:end)-p.hopf.t(1:end-1); h=[h h(end)]; 
      for i=2:tl-1; y0d(1:n,i)=(y(1:n,i+1)-y(1:n,i-1))./(2*h(i)); end
      %  y0d(1:n,1)=(y(1:n,2)-y(1:n,end))./h(1);
      y0d(1:n,1)=(y(1:n,2)-y(1:n,1))./h(1); % just forward
      y0d(1:n,end)=0*y0d(1:n,1); %y0d(:,end) not used
      
end
