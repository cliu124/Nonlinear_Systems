function out=cmcbra(p,u)
X=p.X; n=cross(p.mat.Dx*X,p.mat.Dy*X,2);
A=sum(p.mat.M*(sqrt(dot(n,n,2))));
V=sum(p.mat.M*(dot(X,n,2)))/3;
if V<0; V=-V; end
out=[u(p.nu+1:end);A;V]; 