function p=oosetfemops(p) % dx^2 via cheb 
n=p.np; p.mat.M=speye(n); g=[0 1 0; 0 1 0]; % NBCs 
[xt,D2]=cheb2bc(n,g); % following Reddy-Weideman 
p.mat.L=D2/(p.lx^2); p.x=p.lx*xt; % rescaled Laplacian and mesh 