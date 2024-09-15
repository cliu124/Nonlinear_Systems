function o=asks(s,i); 
% ASKS: little function to ask user for string, with default 
%
%  o=asks(s,i)
% s=string, i=default
as=[s,'(default=',i,'):'];
reply=input(as, 's'); 
if isempty(reply) reply=i; end
o=reply;
end
