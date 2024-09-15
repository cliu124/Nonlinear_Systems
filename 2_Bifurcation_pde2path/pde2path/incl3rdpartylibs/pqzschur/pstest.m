%% a very simple pqzschur-test: choose n,m 
n=4; m=2; AA=zeros(n,n,m);
for i=1:m; AA(:,:,i)=rand(n)+1i*rand(n); end 
BB=AA;
%% with return of Q,Z 
[A,B,Q,Z] = pqzschur(AA,BB);
%% without return of Q,Z 
[A,B] = pqzschur(AA,BB);
