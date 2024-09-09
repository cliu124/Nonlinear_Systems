function Gu=getGupde(p,u,r)
% GETGUPDE: get jacobian for pde part used in full jacobian of getGu.
%
%  Gu=getGupde(p,u,r)
% Here "r" is the residual at u used in numjac
% Imporant switches: p.sw.jac, p.sw.sfem
%
% See also getGu, filltrafo, stanparam, numjac
global pj; % numjacfac numjacG; % for call of numjac
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
else  thresh=p.nc.del*ones(p.nu,1); pj=p; pj.u=u; 
    M=getM(p); 
    try spmat=p.mat.spmat; 
    catch spmat=ones(p.nc.neq); %fprintf('p.mat.spmat not set, using block all-to-all\n'); 
    end
    np=p.nu/p.nc.neq; M=kron(spmat,M(1:np,1:np));  S=ones(p.np); 
    u=u(1:p.nu); %spy(S); pause 
    [Gu,njfac,njG]=numjac('resinj',0,u(:),r(1:p.nu),thresh,[],0,S,[]); 
end