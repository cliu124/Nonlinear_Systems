function c=isoc(divmat,neq,nt)
% ISOC: c-tensor for the isotropic case
%
% * cij12=cij21=0 (no (x,y)-mixed derivative terms) 
% * cij11=cij22 (isotropic diffusion)
%
%  c=isoc(divmat,neq,nt)
%
% Hence, c is only a N-by-N tensor, where the i-th component of 
% div(c tensor grad u) is div(divmat_i dot grad u),
% where divmat_i ist i-th row of divmat.
%  - divmat is a neq-by-neq matrix of vectors with nt components
%  - In particular, if diffusion/fluxes are space-independent then use
%    nt=1 and divmat an neq-by-neq matrix.
%
% See also assempde.  Legacy setup, somewhat obsolete. 
c=[];
dz=zeros(1,nt);
for j=0:neq-1
    for i=1:neq
        cij11=divmat(i,1+j*nt:(j+1)*nt);
        cij22=cij11;
        cij12=dz;
        cij21=dz;
        c=[c;cij11;cij21;cij12;cij22];
    end
end

% ordering example: (3 compo system, i.e. 4*9=36 entries in c) 
%c=[c1111; c1121; c1112; c1122; c2111; c2121; c2112; c2122; ...
%   c3111; c3121; c3112; c3122; ...
%   c1211; c1221; c1212; c1222; c2211; c2221; c2212; c2222;...
%   c3211; c3221; c3212; c3222; ...
%   c1311; c1321; c1312; c1322; c2311; c2321; c2312; c2322;...
%   c3311; c3321; c3312; c3322];
