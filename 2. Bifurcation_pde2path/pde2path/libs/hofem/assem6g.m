function [K,M]=assem6g(p,varargin)
% assem6g: assemble K and M for 6nodes triangles, general coefficients 
% Adapted from Pozrikidis' FSElib for x-dependent and anisotr. diffusion 
t=p.hofem.tri; po=p.pdeo.grid.p'; ng=size(po,1); gdm=0*speye(ng); gmm=gdm;  
a=1; c=1; isw=1; try isw=p.hofem.isw; catch; end % switch for interpol of c 
switch nargin; 
    case 2; c=varargin{1}; % diffusion tensor 
    case 3; c=varargin{1}; a=varargin{2}; % nodal values of diff and mass
end 
if size(a,1)==1; a=a*ones(ng,1); end 
[n,m]=size(c); % n,m, pause 
if n==1; 
   if m==1; c1=c(1,1)*ones(1,ng); c2=c1;  % just a scalar 
   else % c(x)*(\pa_x^2+\pa_y^2) 
       c1=c(1,1:ng); c2=c1; 
   end
else % c is matrix
    if m==2; ong=ones(1,ng); c1=c(1,1)*ong; c4=c(2,2)*ong; c2=c(1,2)*ong; c3=c(2,1)*ong;% const coeff
    else c1=c(1,1:ng); c2=c(1,ng+1:2*ng); c3=c(2,1:ng); c4=c(2,ng+1:2*ng); % matrix function 
    end
end 
ne=size(t,1); NQ=12; try NQ=p.hofem.nq; end; 
for l=1:ne     % loop over the elements
% compute the element diffusion and mass matrices
j=t(l,1); x1=po(j,1); y1=po(j,2); ce1(1)=c1(j); ce2(1)=c2(j); ce3(1)=c3(j); ce4(1)=c4(j); 
j=t(l,2); x2=po(j,1); y2=po(j,2); ce1(2)=c1(j); ce2(2)=c2(j); ce3(2)=c3(j); ce4(2)=c4(j); 
j=t(l,3); x3=po(j,1); y3=po(j,2); ce1(3)=c1(j); ce2(3)=c2(j); ce3(3)=c3(j); ce4(3)=c4(j); 
j=t(l,4); x4=po(j,1); y4=po(j,2); ce1(4)=c1(j); ce2(4)=c2(j); ce3(4)=c3(j); ce4(4)=c4(j); 
j=t(l,5); x5=po(j,1); y5=po(j,2); ce1(5)=c1(j); ce2(5)=c2(j); ce3(5)=c3(j); ce4(5)=c4(j); 
j=t(l,6); x6=po(j,1); y6=po(j,2); ce1(6)=c1(j); ce2(6)=c2(j); ce3(6)=c3(j); ce4(6)=c4(j); 
ce=[ce1 ce2; ce3 ce4]; 
[edm_elm, emm_elm, arel]=edmm6g(x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, NQ,ce,isw);

  for i=1:6 
     i1=t(l,i); % row of gdm 
     for j=1:6
       j1=t(l,j); 
      gdm(i1,j1)=gdm(i1,j1) +edm_elm(i,j);
      gmm(i1,j1)=gmm(i1,j1) +emm_elm(i,j);
     end
   end
end
K=sparse(gdm); M=sparse(gmm); 
