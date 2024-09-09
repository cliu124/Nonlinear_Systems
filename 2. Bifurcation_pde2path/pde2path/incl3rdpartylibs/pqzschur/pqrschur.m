function [Q,A] = pqrschur(A);

Aold = A;

p = size(A,3);
n = size(A,1);
B = zeros(n,n,p);
for ind = 1:p,
   B(:,:,p-ind+1) = A(:,:,ind);
end
s = ones(1,p);
if isreal(B),
   B = complex(B);
end
[QX,B] = percomplex(2,B,s);
Q = zeros(n,n,p);
for ind = 1:p,
   A(:,:,p-ind+1) = B(:,(ind-1)*n+1:ind*n);
end
for ind = 2:p,
   Q(:,:,p-ind+2) = QX(:,(ind-1)*n+1:ind*n);
end
Q(:,:,1) = QX(:,1:n);

err = norm(Q(:,:,1)'*Aold(:,:,p)*Q(:,:,p) - A(:,:,p)) / norm(A(:,:,p));
err = max(err,norm(A(:,:,p) - triu(A(:,:,p))) / norm(A(:,:,p)));
for ind = 2:p,
   err = max(err,norm(Q(:,:,ind)'*Aold(:,:,ind-1)*Q(:,:,ind-1) - A(:,:,ind-1)) / norm(A(:,:,ind-1)));
   err = max(err,norm(A(:,:,ind-1) - triu(A(:,:,ind-1))) / norm(A(:,:,ind-1)));
end
if err > 10^(-14),
   error('pqrschur failed!');
end
