function p=initeig(p,varargin)
% INITEIG: init eigref=guesses for imag parts of evals 
%
%  p=initeig(p)  uses estimate ommax=||G_u|| for cutoff on the imag axis 
%  p=initeig(p, ommax)   uses user-supplied ommax 
% 
if ~isfield(p,'hopf'); p=initwn(p,1,2); end  % init test vectors Schur complement
r=resi(p,p.u); Gu=getGu(p,p.u,r);
if nargin>1; ommax=varargin{1}; else; ommax=norm(Gu,1); end 
kv=[0 1/32 1/16 1/10 1/8 3/16 1/4 3/8 1/2 3/4]; kv=ommax*kv; kl=length(kv); 
nj=p.hopf.wn.j; g=zeros(nj,5*kl); n=size(Gu,1);  
if(p.hopf.wn.Msw==0) un=speye(n,n); 
else un=p.mat.M; end 
if(p.nc.nq>0) % trivially extend M in case of auxiliary equations
    un=[[un zeros(p.nu,p.nc.nq)]; [zeros(p.nc.nq,p.nu) speye(p.nc.nq)]]; 
end
for j=1:nj % loop over test vectors 
  b=p.hopf.wn.b(:,j); c=p.hopf.wn.c(:,j); 
  go=c'*(Gu\b); gov=[real(go);imag(go)]; g(j,1)=go;
  l=2; jv(j)=0; 
  while l<kl;  ok=0; 
    while ~ok; z=kv(l)*1i; 
        gn=c'*((Gu'-z*un)\b); 
        gnv=[real(gn);imag(gn)]; 
      ca=gov'*gnv/(norm(gov)*norm(gnv)); 
      if ca<0.7 % angle too large, refine kv 
        kref=[(kv(l-1)+kv(l))/2, kv(l), (kv(l)+kv(l+1))/2]; 
        kv=[kv(1:l-1), kref, kv(l+1:kl)]; kl=kl+2;
      else ok=1; end
    end
    g(j,l)=gn; gov=gnv; l=l+1; 
  end
end
if nj>1; ga=max(abs(g)); else ga=abs(g); end; 
[ga,idx]=sort(ga,'descend'); omh=kv(idx(1)); p.hopf.wn.gh=ga(1);
%figure(7); clf; plot(real(g(j,:)),imag(g(j,:)),'*-'); pause; 
figure(8); clf; plot(kv(1:kl-1),abs(g(j,1:kl-1)),'-*'); 
set(gca,'FontSize',14); 
if omh==0; fprintf('no om-Hopf found, doing nothing.\n'); % p.nc.nom=1; p.nc.omv=p.nc.eigref; 
else; fprintf('setting om(2)=%g\n', omh); 
p.nc.eigref=[0, omh*1i]; p.nc.neig=[4, 4]; % p.sol.restart=1; % restart
if p.sol.restart~=1; % (re)init of ineg(j), j>=2
  for j=1:length(p.nc.eigref); [p.sol.ineg(j),muv,V]=spcalc(Gu,p,j); end; end
end