function y=myeigsfu(p,A,B,sig,b) % solver to be called in eigs; calls lssgmres 
% (with global prec); here for SH: typical matrix free spectral differentiation 
global p2pglob; fu=p2pglob.fu; p2pglob.fu=fu+sig; % put shift into fu
y=lssgmres(A,b,p); p2pglob.fu=fu;  % solve and restore fu 