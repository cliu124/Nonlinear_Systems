function mymov2avi(mov,fn)
% mymov2avi: HU-interface to VideoWriter (which seems strange under linux!) 
v=VideoWriter(fn,'Archival'); v.FrameRate=3; open(v); 
for i=1:length(mov); writeVideo(v,mov(i)); end; close(v); 
fprintf('run in shell :\n'); 
%str=['avconv -i ' fn '.mj2 -q 10 ' fn '.avi; rm ' fn '.mj2; vlc ' fn '.avi']; disp(str); 
str=['ffmpeg -i ' fn '.mj2 -q 4 ' fn '.avi; vlc ' fn '.avi']; disp(str); 
fprintf('or use some other player.\n'); 
end 
%ffmpeg -i m1.avi.mj2 m1.avi
