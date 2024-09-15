function r=sGpf(p,u) % for brusselator with time periodic forcing  
f=nodalfpf(p,u); par=u(p.nu+1:end); % a,b,d_u, d_v,al,om
du=par(3); dv=par(4); del=par(5); om=2*pi*par(6); 
K=kron([[du,0];[0,dv]],p.mat.K); v1=u(p.nu-1); v2=u(p.nu); va=v1^2+v2^2; 
ru=K*u(1:p.nu-2)-p.mat.M0*f; 
rv1=-del*v1+om*v2+va*v1; rv2=-om*v1-del*v2+va*v2; % the osc.eqns
r=[ru;rv1;rv2]; % append osc.eqns to standard PDE residual 