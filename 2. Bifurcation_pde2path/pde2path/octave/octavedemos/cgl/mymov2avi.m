function mymov2avi(mov,fn)
% HU mov2avi (VideoWriter seems strange under linux!) 
v=VideoWriter(fn,'Archival'); v.FrameRate=3; open(v); 
for i=1:length(mov); writeVideo(v,mov(i)); end; close(v); 
fprintf('run in shell :\n'); 
str=['avconv -i ' fn '.mj2 -q 10 ' fn '.avi; rm ' fn '.mj2; xplayer ' fn '.avi']; disp(str); 
fprintf('or use some other player.\n'); 
end 
