function p=oosetfemops(p) % ac2D 
gr=p.pdeo.grid; np=p.np; p.sf=1e3; % stiff spring constant for DBC via Robin-BC
[K,M,~]=p.pdeo.fem.assema(gr,1,1,1);  % indep. of BC 
p.mat.K=[K 0*K; 0*K K]; 
switch p.nc.neq
    case 2; p.mat.M=[M 0*M; 0*M M]; 
    case 4; p.mat.M=[M 0*M 0*M 0*M; 0*M M 0*M 0*M;  0*M 0*M M 0*M; 0*M 0*M 0*M M];
    otherwise; fprintf('wrong neq\n'); return;
end
bc1=gr.robinBC(1,0); gr.makeBoundaryMatrix(bc1); 
[p.mat.Q,p.mat.G,~,~]=p.pdeo.fem.assemb(gr); % the BC matrices 
po=getpte(p); x=po(1,:); y=po(2,:); p.mat.pot=pot(p); 
Kr=convection(p.pdeo.fem,p.pdeo.grid,[-y;x]); % convection: 
p.mat.Krot=Kr; %p.mat.Krv=[Kr 0*Kr; 0*Kr Kr]; 
p.mat.R=[0*speye(np), -speye(np); speye(np) 0*speye(np)];  % to freeze gauge 