function p=initwn(p,j,sw)
% initwn: set b vectors to prepare initeig (winding nr comp)
n=p.nu+p.nc.nq; np=p.np; nv=1:n; nv=nv'; 
b1=ones(n,1); b2=rand(n,1); b3=rand(n,1); b4=rand(n,1); b5=rand(n,1);  
switch sw;
    case 1; fprintf('random\n'); 
    case 2; b2=-n/2+nv; b3=(nv-n/2).^3; b4=nv.^2; b5=rand(n,1); % seems best!
    case 3; [po]=getpte(p); x=po(1,:)'; x1=min(x); x2=max(x); lx=x2-x1; 
            b1=ones(np,1); b2=cos(0.5*pi*(x-x1)/lx); b3=cos(pi*(x-x1)/lx); 
            b4=cos(1.5*pi*(x-x1)/lx); b5=cos(2*pi*(x-x1)/lx); size(b5)
            b2=[b1;b2]; b3=[b1;b3]; b4=[b1;b4]; b5=[b1;b5]; b1=[b1;b1]; 
    case 4; a=ones(n/2,1); b1=[a; a];  b2=[2*a; a];  b3=[3*a; a];  b4=[4*a; a];  b5=[5*a; a]; 
end
b=[b1, b2, b3, b4, b5]; b=b(:,1:j); wn.j=j; 
for j=1:wn.j; b(:,j)=b(:,j)/norm(b(:,j)); end 
wn.b=b; wn.c=b; wn.Msw=1; p.hopf.wn=wn; 