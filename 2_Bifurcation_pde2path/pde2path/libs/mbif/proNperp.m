function y=proNperp(x,psiv,phiv)
% proNper: project x onto orthog. complement of span(phi(..))
m=size(psiv,2); xc=0*x; for i=1:m; xc=xc+(x'*psiv(:,i))*phiv(:,i); end; 
y=x-xc; 