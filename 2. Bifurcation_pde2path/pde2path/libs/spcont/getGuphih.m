function Guuphi=getGuphih(p,u,phir,phii)
% getGuTpsiu: \pa_u(G_u'*\psi) via numjac 
% r=Gu'*psi 
global pGu; 
try Xcont=p.sw.Xcont; catch Xcont=0; end 
try njt=p.nc.njthreshsp; catch njt=1e3; end; 
thresh=njt*ones(p.nu,1);  
pGu=p; pGu.u=u; pGu.phir=phir; pGu.phii=phii;  % to pass pars and phi to resiGupsi
if Xcont>0; try fM=p.mat.fill; catch; fM=1; end 
    M=massmatrix(p.X,p.tri,'full'); M=fM'*M*fM;     
else;  M=getM(p); end
r0=pderesi(p,u); Gu=getGupde(p,u,r0); r1=[Gu*phir;Gu*phii]; % reference residual 
%10, psipde(1:3)', u(p.nu+1:end)', u(1:4)', r0(1:4)', r(1:4)', pause 
try spmat=p.mat.spmat; catch spmat=ones(p.nc.neq); end
np=p.nu/p.nc.neq; M=kron(spmat,M(1:np,1:np));  
S=(M~=0); S=[M;M]; v=u(1:p.nu); %fig(12); spy(S), pause    
[Guuphi,njfac,njG,nf1,nf2]=numjac('resiGuphih',0,v,r1,thresh,[],0,S,[]); 
if p.sw.verb>2; %fprintf('numjac used %i calls to ***resiGuphih***\n',nf1); 
        del=numjacinfo('resiGuphih',0,v,r1,thresh,nf1); 
end