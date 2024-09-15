function [K,M] = ass1dlob (ne,xe,np,ng,c)
% assem1dlob: assemble K,M for lobatto-FEM, following FSELIB 
% ne=#elems, xe=end-points, np(1:xe)=polynomial degrees, ng=#global nodes 
% c=connectivity matrix
h=diff(xe); K=zeros(ng,ng); M=K; 
for l=1:ne  % loop over the elements
   m=np(l);
    elm_mm=0.5*h(l)*emm_lob_tbl(m);         % element mass matrix
%   elm_mm=0.5*h(l)*emm_lob_lump_tbl(m);    % element mass matrix
    elm_dm=2.0*edm_lob_tbl(m)/h(l);         % element diffusion matrix
   for ip=1:m+1
      i1=c(l,ip);
      for jp=1:m+1
       i2=c(l,jp);
       K(i1,i2)=K(i1,i2) + elm_dm(ip,jp);
       M(i1,i2)=M(i1,i2) + elm_mm(ip,jp);
      end
   end
end 