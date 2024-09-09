function p=nlbinit(p,lx,nx,om,kstar,sw_asymp)
p.file.dircheck=0; dir=sprintf('%s',inputname(1)); p=setfn(p,dir);
p.fuha.G=@nlbf; p.fuha.Gjac=@nlbjac; p.fuha.postmmod=@nlbpmm; % only needed for mesh-ref. 
p.fuha.sG=@nlbsG; p.fuha.sGjac=@nlbsjac; p.fuha.outfu=@nlbbra;
ly=lx; [p.mesh.geo,bc]=recnbc2(lx,ly);  %zero Neumann BC (later replaced by periodic BC)
p.fuha.bc=@(p,u,lam) bc; p.fuha.bcjac=@(p,u) bc;
ny=nx; p=stanmesh(p,nx,ny); p=setbmesh(p); p.sol.xi=1/(2*p.nu); 
p.mat.pot=pot(p); p.mat.poti=pdeintrp(p.mesh.p,p.mesh.t,p.mat.pot); %periodic potential
p.u=zeros(p.nu+3,1);

p=box2per(p,[1,2]); % generate matrices for periodic BCs, reduce sol. vector

p.u=zeros(p.nu+3,1); p.nc.ilam=1; % p.nu, size(p.u)
p.eqn.c=[1;0;0;1;1;0;0;1]; b=zeros(8,1); b(5)=-2*kstar(1);
b(3)=2*kstar(1);b(6)=-2*kstar(2); b(4)=2*kstar(2); p.eqn.b=b;
p.u(p.nu+1)=om; p.u(p.nu+2)=kstar(1); p.u(p.nu+3)=kstar(2);
p.eqn.a=0; p.sw.sfem=0; %p=setfemops(p);
%-----------------------------------------------------
%initial guess based on asymptotics
if sw_asymp==1
    %compute the linear Bloch wave at k =kstar near a chosen omega=om point
    p.u(p.nu+1)=0;   %set omega=0  for the linear eigenvalue problem
    r=pderesi(p,p.u); Gu=getGupde(p,p.u,r);
    p.nc.eigref=om; [ineg,muv,V]=spcalc(Gu,p);
    fprintf('first omegas:\n');
    muv(1:8), fprintf('pick one (real) by index\n');
    evnr=1; evnr=asknu('nr (0 for end)',evnr); if evnr==0; return; end
    Phi=V(:,evnr); omega= muv(evnr);
    u_full=p.mat.fill*Phi;
    c1=sqrt(triint((u_full(1:p.np).^2+u_full(p.np+1:2*p.np).^2),...
        p.mesh.p,p.mesh.t));
    p.ev=Phi/c1;
    c2=triint((u_full(1:p.np).^2+u_full(p.np+1:2*p.np).^2).^2,...
        p.mesh.p,p.mesh.t);
    %amplitude in the asymptotic approximation of the solution
    A=sqrt(1/c2); eps=0.1;  %too small eps results in convergence to 0
    p.u(1:p.nu)=A*eps*p.ev; plotsol(p,1,1,3);plotsol(p,6,2,3);
    p.u(p.nu+1)=omega+sign(p.sig)*eps^2;
    r=norm(pderesi(p,p.u),2); fprintf('initial res=%g\n',r);
    [p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res);
    plotsol(p,1,1,3);plotsol(p,6,2,3);
end
p.tau=zeros(length(p.u),1); p.tau(length(p.u))=1;
