function cpplot(poc,ncomp,varargin)
ts=4; if nargin>2; ts=varargin{1}; end 
cp=poc.cp; v=[70,30]; s1=poc.oc.s1; 
x=getpte(s1); x=x(1,:); np=s1.np; T=real(cp.par); tv=T*cp.t(1:ts:end); 
u=cp.u; ga=s1.u(s1.nu+5); 
for k=1:ncomp+1 % x-t-plots of solns   
    figure(k); clf; hold on;    
    if k<=ncomp; z=u((k-1)*np+1:k*np,1:ts:end)'; % states and co-states 
    else % control 
        l1=u(2*np+1:3*np,1:ts:end)';
        z=-(1+l1)./ga; % control 
    end
    surf(x,tv,z);
    ylabel('t'); view(v); 
    switch k; case 1; zlabel('v');  case 2; zlabel('w'); 
    case 3; zlabel('\lambda');  case 4; zlabel('\mu'); case 5; zlabel('q'); 
    end 
    set(gca,'fontsize',14); axis tight; 
end