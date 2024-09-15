function [p,r,Gu,Glam,mu,mua]=getker(dir,fname,varargin)
% getker: get kernel of G_u, used in q(c)swibra, or as prep for gentau 
m=1;
if ischar(dir); p=loadp(dir,fname);
else; p=dir; if nargin>1; m=fname; end % input is struct, fname is aux 
end
u=p.u; 
if ~isempty(varargin); m=varargin{1}; end 
r=resi(p,u); [Gu,Glam]=getder(p,u,r);  % residual and jacobian
M=getM(p); %p.mat.M; 
if(p.nc.nq>0) % trivially extend M in case of auxiliary equations
    if isfield(p.fuha,'qMfu'); [qL,qU,qD]=p.fuha.qMfu(p); 
    else [qL,qU,qD]=stanqM(p); 
    end 
    M=getM(p); M=[[M(1:p.nu,1:p.nu) qU]; [qL qD]]; 
end
vs=size(Gu,1); opts.v0=ones(vs,1)/vs; opts.disp=0; 
[phiv,mu]=myeigs(Gu,M,m,0,opts,p); phiv=real(phiv);   % getting EVecs
[psiv,mua]=myeigs(Gu',M,m,0,opts,p); psiv=real(psiv); psiv=hugram2(psiv,phiv); % adjoints 
mud=diag(mu); [muds,idx]=sort(abs(real(mud))); mu=mud(idx)';
muad=diag(mua); [muads,idxa]=sort(abs(real(muad))); muad=muad(idxa)'; 
p.mat.phiv=phiv(:,idx); p.mat.psiv=psiv(:,idx);