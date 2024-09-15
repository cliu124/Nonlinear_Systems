function out=meshqdat(p)
% meshqdat: mesh-quality data; 
% max(h/r), max(A),min(A), max(h), min(h) 
X=p.X; T=p.tri; 
A=doublearea(X,T)/2; e=edge_lengths(X,T);
el=max(e,[],2); % long edges 
s=sum(e,2)/2; r=A./s; q=max(el./r); % inradii.and "quality" 
% R=e(:,1).*e(:,2).*e(:,3)./(8*A);  % outradius
a1=max(A); a2=min(A); 
h1=max(el); h2=min(el); 
out=[q; a1; a2; h1; h2]; 