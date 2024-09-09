function hogradinf(p) 
% hogradinf: print info on |u'|_infty
ho=p.hopf; y0d=ho.y0d; cmax=max(abs(y0d)); t=ho.T*ho.t; tl=ho.tl; 
fprintf('i:    '); for i=1:tl; fprintf('%8i ',i); end; fprintf('\n');
%fprintf('t:        '); for i=1:tl; fprintf('%s ', printcon(t(i))); end; fprintf('\n');
fprintf('t:         '); for i=1:tl; fprintf('%6.2f   ', t(i)); end; fprintf('\n');
%fprintf('t:        '); for i=1:tl; fprintf('%s ', printcon(cmax(i))); end; fprintf('\n');
fprintf('max(udot):'); for i=1:tl; fprintf('%6.2e ', cmax(i)); end; fprintf('\n');
[cmaxs,idx]=sort(cmax,'descend');
fprintf('max of |udot| at :\n'); ni=8; 
fprintf('i:    '); for i=1:ni; fprintf('%8i ',idx(i)); end; fprintf('\n');
%fprintf('t:        '); for i=1:tl; fprintf('%s ', printcon(t(i))); end; fprintf('\n');
fprintf('t:         '); for i=1:ni; fprintf('%6.2f   ', t(idx(i))); end; fprintf('\n');
%fprintf('t:        '); for i=1:tl; fprintf('%s ', printcon(cmax(i))); end; fprintf('\n');
fprintf('max(udot):'); for i=1:ni; fprintf('%6.2e ', cmax(idx(i))); end; fprintf('\n');
end 