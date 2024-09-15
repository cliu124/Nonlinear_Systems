function lfplottip(varargin)
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; lev=varargin{4}; anf=5;
else  wnr=varargin{2}; lev=varagin{3}; anf=4;
end
try aux=varargin{anf}; catch; aux=[]; end
try tol=aux.tol; catch; tol=lev(1)/2; end  
s=p.u(p.nu+p.spar); 
y=p.hopf.y; tl=p.hopf.tl; T=p.hopf.T; tv=T*p.hopf.t; om1=s/(2*pi); om2=1/T; 
try pertol=aux.pertol; catch; pertol=1e-1; end 
try mr=aux.mr; catch mr=0; end
if mr==0 % find m s.t. orbit is (approx) mT per in lab-frame 
  for i=1:20; mr=i; mt=i*om1/om2; mr=round(mt); d=mt-mr; if abs(d)<pertol; mr=i; break; end; end;
end
mt=mr*om1/om2; d=mt-round(mt); 
%if i==20; mr=i; end 
fprintf('s=%g, RW-freq om1=%g, mod-freq om2=%g, om1/om2=%g, used m=%g, miss=%g\n',s,om1,om2,om1/om2,mr,d); 
xtv=zeros(2,mr*(tl-1)); 
xv0=gettip(p,y(1:p.nu,1),lev,tol); 
for m=1:mr
for i=1:tl-1
    t=(m-1)*T+tv(i); xv=gettip(p,y(1:p.nu,i),lev,tol); %xv(1)=0.5*xv(1); 
    R=[cos(s*t) -sin(s*t); sin(s*t) cos(s*t)];
    xtv(:,(m-1)*(tl-1)+i)=R*(0*xv0'+xv')+0*xv'; 
end
end
figure(wnr); clf; plot(xtv(1,:),xtv(2,:),'-'); hold on;  axis image
plot(xtv(1,1),xtv(2,1),'r*'); plot(xtv(1,end),xtv(2,end),'b*');
set(gca,'FontSize',p.plot.fs);