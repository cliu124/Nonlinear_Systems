function [f,J]=huassem3(p,y,T,par) % backward FDs
mxsw=0;  % for 1, use M from X=X+u.*N 
try Xsw=p.sw.Xcont; catch; Xsw=0; end
M=getM(p);  idx=find(all(M==0,2)); % algebraic compos
t=p.hopf.t; tl=length(t); h=diff(t); 
nu=p.nu; na=tl*nu; f=zeros(tl*nu,1); 
% i=1: i-1 at i=tl-1
i=1; y1=[y(:,i);par]; 
if Xsw>0; p.i=i; X=p.hopf.X(:,:,i);
  if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end
  M=getM(p,X)/h(i);  
end
G=p.fuha.sG(p,y1); 
G(idx)=G(idx)/T; % lag
f((i-1)*nu+(1:nu))=M*(y(:,i)-y(:,tl-1))+T*G; 
for i=2:tl-1 % main block 
    y1=[y(:,i);par]; 
    if Xsw>0; p.i=i; X=p.hopf.X(:,:,i);
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
    end
    G=p.fuha.sG(p,y1); 
    G(idx)=G(idx)/T; % lag
    f((i-1)*nu+(1:nu))=M*(y(:,i)-y(:,i-1))+T*G; 
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1); % periodicity 
if nargout==1; return; end
J=sparse(na); 
i=1; % 1st block 
si=(i-1)*nu; y1=[y(:,i); par]; 
if Xsw>0; p.i=i; X=p.hopf.X(:,:,i);
  if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
  M=getM(p,X)/h(i); 
end
Gu=getGupde(p,y1); Gu(idx,:)=Gu(idx,:)/T; 
J(si+(1:nu),si+(1:nu))=T*Gu+M; 
J(si+(1:nu),si+(tl-2)*nu+(1:nu))=-M;
for i=2:tl-1 % 
    y1=[y(:,i); par]; si=(i-1)*nu; 
    if Xsw>0; p.i=i; X=p.hopf.X(:,:,i);
       if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
       M=getM(p,X)/h(i); 
    end
    Gu=getGupde(p,y1); Gu(idx,:)=Gu(idx,:)/T; 
    J(si+(1:nu),si+(1:nu))=T*Gu+M;     
    J(si+(1:nu),si-nu+(1:nu))=-M; 
end 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu);
end 