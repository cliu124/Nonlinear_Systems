function [f,J]=huassem4(p,y,T,par) % like TOM 
%
% [F,JAC,info]=huassem(p,T+del,y,par); 
mxsw=0;  % for 1, use M from X=X+u.*N  
try Xsw=p.sw.Xcont; catch; Xsw=0; end
t=p.hopf.t; tl=length(t); tt=t; h=diff(tt); np=p.np-1; 
nu=p.nu; na=tl*nu; f=zeros(tl*nu,1); 
M=getM(p);  idx=find(all(M==0,2)); % algebraic compos
% i=1: i-1 at i=tl-1
i=1; y1=[y(:,i);par]; p.i=i; 
if Xsw>0; X=p.hopf.X(:,:,i);
  if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
  M=getM(p,X)/h(i);  
end
G=p.fuha.sG(p,y1); G(idx)=2*G(idx)/T; % alg compos
y1=[y(:,tl-1);par];  p.i=tl-1; 
Gm=p.fuha.sG(p,y1); Gm(idx)=0; % zero out alg.compos 
f((i-1)*nu+(1:nu))=M*(y(:,i)-y(:,tl-1))+0.5*T*(G+Gm); 
for i=2:tl-1 % main block 
    y1=[y(:,i);par]; p.i=i; 
    if Xsw>0; X=p.hopf.X(:,:,i);
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
      M=getM(p,X)/h(i);  
    end
    G=p.fuha.sG(p,y1);  G(idx)=2*G(idx)/T; 
    y1=[y(:,i-1);par]; p.i=i-1;
    Gm=p.fuha.sG(p,y1); Gm(idx)=0; % zero out algebraic rows 
    f((i-1)*nu+(1:nu))=M*(y(:,i)-y(:,i-1))+0.5*T*(G+Gm); 
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1); % periodicity 
if nargout==1; return; end
J=sparse(na); 
i=1; % 1st block 
si=(i-1)*nu; y1=[y(:,i); par]; p.i=i; 
if Xsw>0; X=p.hopf.X(:,:,i);
   if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
   M=getM(p,X)/h(i); 
end
Gu=getGupde(p,y1); Gu(idx,:)=2*Gu(idx,:)/T; 
y1=[y(:,tl-1);par];  p.i=tl-1; 
Gum=getGupde(p,y1); Gum(idx,:)=0; % zero out alg rows 
J(si+(1:nu),si+(1:nu))=0.5*T*Gu+M; 
J(si+(1:nu),si+(tl-2)*nu+(1:nu))=0.5*T*Gum-M; 
for i=2:tl-1 %     
    y1=[y(:,i); par]; si=(i-1)*nu; p.i=i; 
    if Xsw>0;  X=p.hopf.X(:,:,i);
      if mxsw; N=getN(p,X); X=X+y(1:np,i).*N; end 
      M=getM(p,X)/h(i); 
    end
    Gu=getGupde(p,y1); Gu(idx,:)=2*Gu(idx,:)/T; 
    i=i-1; y1=[y(:,i); par]; p.i=i; 
    Gum=getGupde(p,y1); Gum(idx,:)=0; %Gu(end,:); 
    J(si+(1:nu),si+(1:nu))=0.5*T*Gu+M;     
    J(si+(1:nu),si-nu+(1:nu))=0.5*T*Gum-M; 
end 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu);
end 