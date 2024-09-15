function [ineg,muv,V]=spcalc(Gu,p,varargin) 
% SPCALC: compute p.nc.neig(j) EVals of Gu near p.nc.eigref(j), j=varargin 
%
%  [ineg,muv]=spcalc(Gu,p)
%  [ineg,muv,V]=spcalc(Gu,p) - also compute and return eigenvectors in V
% 
% * muv=  eigenvalues sorted wrt real parts,
% * ineg= number of negative evals.
if (p.sw.spcont~=0)
    [p,spcont,prim]=spreduce(p); 
    Gu=Gu(1:p.nu+p.nc.nq,1:p.nu+p.nc.nq);
end
vs=size(Gu,1); M=getM(p); 
j=1; if nargin>2; j=varargin{1}; end 
om=p.nc.eigref(j); neig=p.nc.neig(j); 
if(p.nc.nq>0) % trivially extend M in case of auxiliary equations
    if isfield(p.fuha,'qMfu'); [qL,qU,qD]=p.fuha.qMfu(p); 
    else [qL,qU,qD]=stanqM(p); 
    end 
    M=[[M(1:p.nu,1:p.nu) qU]; [qL qD]]; 
end
%mclf(23); spy(M); pause 

if(strcmp(p.sw.eigmeth,'eigs')) % compute p.nc.neig FEM eigenvalues 
 opts.v0=ones(vs,1)/vs; opts.disp=0;  [V,mu]=myeigs(Gu,M,neig,om,opts,p);
    muv=diag(mu); 
else % compute Evals with real-part in [p.nc.eigint(1), p.nc.eigint(2)], usually slow 
    evalc('[V,muv]=sptarn(Gu,M,p.nc.eigint(1),p.nc.eigint(2));');
end
mul=length(muv); [ms,ix]=sort(abs(real(muv))); Vs=zeros(vs,mul);
ineg=0; mus=zeros(1,mul); 
try intol=p.nc.intol; catch intol=0; end 
for i=1:mul; mus(i)=muv(ix(i)); Vs(:,i)=V(:,ix(i));
    if (real(mus(i))<intol) ineg=ineg+1; end
end; 
if p.sw.verb>1; figure(5+j); clf; plot(real(mus), imag(mus),'*'); grid on; end 
minr=min(real(muv)); maxr=max(real(muv)); 
if (p.sw.verb>2 && abs(minr)>abs(maxr)/2)
    fprintf('warning: mu-min=%g, mu-max=%g; consider increasing p.nc.neig\n', mus(1),mus(mul));
end
muv=mus; V=Vs; 

if (p.sw.spcont~=0); p=spextend(p,spcont,prim); end % reset to spectral cont
end

