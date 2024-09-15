function [f,J]=huassem(p,T,y,par)
% tomassem: assemble rhs and jac for hopf (arclength), extracted from TOM
%
% [F,JAC,info]=huassem(p,T+del,y,par); 
t=p.hopf.t; tl=length(t); tt=t; h=diff(tt); %h=1/tl*ones(1,tl-1); %h, pause 
nu=p.nu; na=tl*nu; 
f=zeros(tl*nu,1);% nu, tl, size(f), pause 
for i=1:tl-1 % forward FDs 
    y1=[y(:,i);par]; 
    G=p.fuha.sG(p,y1); M=getM(p)/h(i);     
    f((i-1)*nu+(1:nu))=M*(y(:,i+1)-y(:,i))+T*G;  
end 
f((tl-1)*nu+(1:nu))=y(:,tl)-y(:,1);
if nargout==1; return; end
J=sparse(na); %h, pause 
for i=1:tl-2 % forward FDs     
    y1=[y(:,i); par]; si=(i-1)*nu;   
    Gu=getGupde(p,y1); 
    M=getM(p)/h(i); %i, full(M(1:5,1:5)), pause 
    J(si+(1:nu),si+(1:nu))=-M+T*Gu; 
    J(si+(1:nu),si+(nu+1:2*nu))=M; 
end 
si=(tl-2)*nu; i=tl-1; % penultimate row: 
y1=[y(:,i); par]; M=getM(p)/h(i); Gu=getGupde(p,y1);     
J(si+(1:nu),si+(1:nu))=-M+T*Gu; 
J(si+(1:nu),1:nu)=M; 
si=(tl-1)*nu; % last row 
J(si+(1:nu),1:nu)=-speye(nu); J(si+(1:nu),si+(1:nu))=speye(nu);
end 