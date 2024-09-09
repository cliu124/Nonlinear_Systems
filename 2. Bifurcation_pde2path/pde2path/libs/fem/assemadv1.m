function B=assemadv1(p,t,b)
% ASSEMADV1: assemble advection matrix for scalar problem.
%
%  B=assemadv1(p,t,b) 
% b*grad u, p=points, t=triangles 
% u=scalar, hence b=2 x m vector with m=1 or m=nt=#triangles.
%
% See also assemadv
it1=t(1,:);it2=t(2,:);it3=t(3,:);np=size(p,2);
[ar,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t); % triangle-areas and base-fu-gradients 
b1=b(1,:); b2=b(2,:); 
% the local entries; thx to Uwe Prüfert for the formulas
E1=ar.*(b1.*g1x+b2.*g1y)*(1/3); 
E2=ar.*(b1.*g2x+b2.*g2y)*(1/3); 
E3=ar.*(b1.*g3x+b2.*g3y)*(1/3);
B=sparse(it1,it1,E1,np,np);
B=B+sparse(it1,it2,E2,np,np);
B=B+sparse(it1,it3,E3,np,np);
B=B+sparse(it2,it1,E1,np,np);
B=B+sparse(it2,it2,E2,np,np);
B=B+sparse(it2,it3,E3,np,np);
B=B+sparse(it3,it1,E1,np,np);
B=B+sparse(it3,it2,E2,np,np);
B=B+sparse(it3,it3,E3,np,np);
