function q=qrot(p,u) % rotational phase condition 
n=p.np; q=(p.mat.Krot(1:n,1:n)*p.u(1:n))'*u(1:n); 