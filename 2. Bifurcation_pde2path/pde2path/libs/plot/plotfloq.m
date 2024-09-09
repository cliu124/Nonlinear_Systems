function h=plotfloq(pd,varargin)
% plotfloq: plot Floquet multipliers in p (p2p-struct, or directory) 
noa=nargin-1; k=1; 
if(ischar(pd)) 
    if (noa>0 && ischar(varargin{1})); p=loadp(pd,varargin{1}); k=k+1; noa=noa-1; 
    else p=loadp(pd); 
    end
else p=pd; 
end
wnr=7; if noa>0; wnr=varargin{k}; end
muv1=p.hopf.muv1;muv2=p.hopf.muv2;
figure(wnr); clf; set(gca,'FontSize',12); 
h=plot(real(muv1), imag(muv1),'*b', real(muv2), imag(muv2), '*r'); 
set(h, 'Markersize',8);
hold on;  s=0:0.1:2*pi+0.1; plot(cos(s),sin(s)); axis image
lam=getlam(p); 
fprintf('lam=%g, ind(u_H)=%i\n',lam,length(muv2)); 

 