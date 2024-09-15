function s=printcon(x)
% PRINTCON: give 10char to number, default %f
%
%  s=printcon(x)
%
% (used to give unified and compact formatting for printout) 
s=num2str(x,'%+9.5f');si=size(s,2);
if (si>10 || abs(x)<1e-5) s=num2str(x,'%+2.3e'); end
if(s(1)=='+') s(1)=' ';end
si=size(s,2);
if si<10; sf=''; % align left
  for i=1:10-si; sf=[sf ' ']; end; 
  s=[sf s]; end
end