function n=hunorm(u)
% hunorm: l2-norm, dim-independent
u2=u.*u; u2=u2(:); n=sqrt(sum(u2));