%% make movie of last points on h2-h6:  
p2=loadp('rw2','pt10'); p3=loadp('rw3','pt10'); 
p5=loadp('rw5','pt10'); p6=loadp('rw6','pt10'); p7=loadp('rw7','pt10');
%%
cnr=1; wnr=4; mov=homovplot5(p2,p3,p5,p6,p7,cnr,wnr); 
%movie2avi(mov,'m5.avi');
%% make movie of 'spirals' (r=3) 
p2=loadp('s2','pt29'); p3=loadp('s3','pt34'); p4=loadp('s5','pt20'); 
p5=loadp('s6','pt18'); p6=loadp('s7','pt18');
%%
mov=homovplot5(p2,p3,p4,p5,p6,cnr,wnr); 
%movie2avi(mov,'m5s.avi');
% on linux: convert with mencoder in.avi -o out.avi -speed 0.4 -ovc lavc 
% or (for web) avconv -i m1.avi -q 10 -vf 'setpts=8*PTS' m1.ogv 
% or    ffmpeg -i in.avi -vf 'setpts=6*PTS' out.avi
% to about 1/100 in size
%% rotating rose! 
p=loadp('s2d','pt44'); mov=homovplot(p,1,1,3); movie2avi(mov,'rose.avi'); 
