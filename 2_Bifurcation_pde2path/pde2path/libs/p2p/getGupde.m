function Gu=getGupde(p,u,r)
% GETGUPDE: get jacobian for pde part used in full jacobian of getGu.
%
%  Gu=getGupde(p,u,r)
% Here "r" is the residual at u used in numjac
% Imporant switches: p.sw.jac, p.sw.sfem
%
% See also getGu, filltrafo, stanparam, numjac
global pj; % for call of numjac
try Xcont=p.sw.Xcont; catch Xcont=0; end 
if (p.sw.jac>0); % pde-part anal. 
    if p.sw.sfem==0 % use full jac
    bc=p.fuha.bcjac(p,u); upde=p.mat.fill*u(1:p.nu); 
    [cj,aj,bj]=p.fuha.Gjac(p,u); zerov=zeros(p.nc.neq,1); 
    [Gu,dum]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,cj,aj,zerov, upde); 
    if(any(bj)); Kadv=assemadv(p.mesh.p,p.mesh.t,bj); 
        Gu=Gu-Kadv; end 
    Gu=filltrafo(p,Gu);
    return;
    else % use "simple jac". Note that p.mat.fill is built in here already
    Gu=p.fuha.sGjac(p,u); 
    end
else pj=p; pj.u=u; % to pass data to resinj     
    try; njt=p.nc.njthresh; catch; njt=p.nc.del; end 
    thresh=njt*ones(p.nu,1); 
    if Xcont>0; try fM=p.mat.fill; catch; fM=1; end 
       M=massmatrix(p.X,p.tri,'full'); M=fM'*M*fM;     
    else;  M=getM(p); end
    try spmat=p.mat.spmat; catch spmat=ones(p.nc.neq); end    
    np=p.nu/p.nc.neq; M=kron(spmat,M(1:np,1:np));  S=(M~=0); v=u(1:p.nu); 
   % fig(12); spy(S), pause   
    [Gu,njfac,njG,nf1,nf2]=numjac('resinj',0,v,r(1:p.nu),thresh,[],0,S,[]);
    if p.sw.verb>3; %fprintf('numjac used %i calls to ***resinj***\n',nf1); 
        del=numjacinfo('resinj',0,v,r(1:p.nu),thresh,nf1); 
    end
   %  fig(12); spy(Gu), pause   
end