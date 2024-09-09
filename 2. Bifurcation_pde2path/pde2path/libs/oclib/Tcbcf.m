function bc=Tcbcf(ua,ub,par,opt)
% Tcbcf: boundary conditions, free T case 
n=length(ua); ugam=opt.u1(1:n);
% regular bc
bc=[ua(1:n/2)-opt.u0(1:n/2);opt.F*(ub-ugam)];

% bc for hopf resp. steady end point
if opt.s1ho==1 % CP to CPS
    if opt.freeT % depreciated, opt.freeT=0!
        bc=[bc;opt.nf*(norm(ub-ugam,2)^2/n-opt.tadev2^2)];
    else
        bc=[bc(1:n/2);opt.F2*(ub-ugam)];
    end
else % CP to CSS
    if opt.freeT
        bc=[bc;(norm(ub-ugam,2)^2/n-opt.tadev2^2)];
        %opt.tadev2, bc(end), pause 
    else
        bc=[bc;par(1)-opt.T];
    end
end
end