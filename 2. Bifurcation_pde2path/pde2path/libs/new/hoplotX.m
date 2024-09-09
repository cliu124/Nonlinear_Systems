function hoplotX(varargin)
[p,anf]=getp(varargin{:}); 
ho=p.hopf; tl=length(ho.t); idx=1:tl; n=p.np;  
noa=nargin; k=anf; fnr=1; lw=6; vi=[50,30]; cmp=2; 
while k<=noa % parse all optional inputs
    ak=lower(varargin{k}); 
    switch ak 
        case 'fnr'; fnr=varargin{k+1}; k=k+2; 
        case 'lw'; lw=varargin{k+1}; k=k+2;
        case 'vi'; vi=varargin{k+1}; k=k+2; 
         case 'cmp'; cmp=varargin{k+1}; k=k+2;    
        otherwise k=k+1; 
    end
end 
y=ho.y; ov=ones(n,1); mclf(fnr); i1=(cmp-1)*n+1; i2=cmp*n; 
for i=idx; 
 X=ho.X(:,:,i); t=ho.T*ho.t(i)*ov; 
% plot3(X(:,1),t,X(:,2),'b'); hold on; 
 patch([X(:,1); nan],[t; nan], [X(:,2); nan],[y(i1:i2,i); nan],... 
    'FaceColor','none','EdgeColor','interp','LineWidth',lw); 
hold on; 
end
ylabel('t');  view(vi); box on; grid on; %xlabel('x'); zlabel('y');
title([p.file.pname mat2str(p.file.count)]); 
set(gca,'fontsize',p.plot.fs); axis tight; 