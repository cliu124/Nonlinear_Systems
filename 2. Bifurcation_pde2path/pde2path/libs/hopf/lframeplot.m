function lframeplot(varargin)
% lframeplot: x-t-plot of rel. Hopf orbit with period T in lab frame 
%
%  lframeplot(p,wnr,cmp,aux) 
%  lframeplot(dir,pt,wnr,cmp,aux) 
% 
% where aux may include: 
%  spar: parameter-index for speed s
%  mper: user provided m s.t. sol. is T^*=mT periodic in lab-frame. 
%  pertol:  tolerance for computation of minimal m
if ischar(varargin{1})
   dir=varargin{1}; pt=varargin{2}; str=[dir,'/',pt,'.mat']; 
   p=loadp(dir,pt); wnr=varargin{3}; cmp=varargin{4}; anf=5;
else  p=varargin{1}; wnr=varargin{2}; cmp=varargin{3}; anf=4;
end
try aux=varargin{anf}; catch; aux=[]; end 
try mper=aux.mper; catch; mper=0; end 
try spar=aux.spar; catch; try; spar=p.spar; catch; spar=0; end; end
try vi=aux.vi; catch; vi=p.plot.view; end 
try cb=aux.cb; catch; cb=0; end 
try pertol=aux.pertol; catch; pertol=1e-1; end 
try s=p.u(p.nu+spar); catch; s=0; end 
z=p.hopf.y; z0=z; tv=p.hopf.t; T=p.hopf.T; tl=length(tv); 
po=getpte(p); x=po(1,1:end-1)'; dx=x(2)-x(1); np=p.np-1; lx=po(1,end)-x(1); 
if mper==0 
  for q=1:20 % find m s.t. orbit is mT per in lab-frame 
    mt=q*lx/(s*T); mr=round(mt); d=mt-mr;   fprintf('d=%g, m=%g\n',d, mr);
    if abs(d)<pertol; mper=mr; break; end
  end; 
end 
tt=[]; for m=1:mper; tt=[tt (m-1)*T+T*tv(1:tl-1)]; end  
fprintf('used m=%g, T^*=%g, last shifts:\n', mper, mper*T); 
for m=0:mper-1
for j=1:tl-1
    j0=m*(tl-1)+j; t=tt(j0); ks=round(s*t/dx); uj=z0((cmp-1)*np+1:cmp*np,j); 
    z((cmp-1)*np+1:cmp*np,j0)=circshift(uj,ks,1); 
end
fprintf('%i, ',rem(ks,np)); 
end
fprintf('\n'); 
sol.x=tt; sol.y=z; 
xtplot(p,sol,wnr,cmp,[10 80],[]); shading interp; figure(wnr); axis image; view(vi); 
if cb; colorbar; end 
ylabel('t'); 
tit=[dir '/pt' mat2str(p.file.count-1)]; %tit=str; 
if anf==5; title(tit); end 

    
    