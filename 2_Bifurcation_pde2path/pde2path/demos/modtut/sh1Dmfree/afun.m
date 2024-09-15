function y=afun(uf) % lin.of SH used in lssgmres; p2pglob.fu filled in sGjac 
global p2pglob; mu=p2pglob.mu; fu=p2pglob.fu;  y=mu.*uf-dct(fu.*idct(uf)); 
