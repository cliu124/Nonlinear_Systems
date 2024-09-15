function y=afun(uf) % A*uf function for lssgmres for SH via fu (from sGjac) 
global p2pglob; mu=p2pglob.mu; fu=p2pglob.fu; F=p2pglob.F; 
y=mu.*uf-F*(fu.*(F'*uf)); 
