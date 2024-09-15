function [z,p]=lssbel(A,w,p)
% LSSBEL: bordered elimination lss after govaerts with post iterations and
% border width p.bel.bw
% 
%  [z,p]=lssbel(A,w,p)
%
% see also: lss, blss, blssbss, bel, belpi, bss
if p.bel.bw>0
 if isequal(p.fuha.innerlss,@lss) % use special bel if inner lss is \ 
     [z,p]=mbellu(A,w,p); return; 
 end 
 n=size(w,1); nm=n-p.bel.bw; 
 [x,y,p,r]=belpi(A(1:nm,1:nm), A(1:nm, nm+1:n), A(nm+1:n, 1:nm), ...
     A(nm+1:n, nm+1:n),w(1:nm),w(nm+1:n),p);
 if r<p.bel.tol; z=[x;y]; 
 else fprintf('bordered elim failed, using innerlss on full system\n'); 
   [z,p]=p.fuha.innerlss(A,w,p);
 end 
else    
  [z,p]=p.fuha.innerlss(A,w,p);  
end
end
