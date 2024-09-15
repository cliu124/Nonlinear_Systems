clc; clear all;%32 is space,10 is newline, 39 is ', 37 is %, 59 is ;
files = dir;
total = 0;
totalI = 0;
for i=3:numel(files)
	if numel(files(i).name) < 2 || not(strcmp(files(i).name(end-1:end),'.m')) || exist(files(i).name) == 7
		continue;
	end;
	fid = fopen(files(i).name);
	A = fread(fid);
	fclose(fid);
	ttl = sum(A == ';'); total = total + ttl; totalI = totalI + 1;
	disp(sprintf('%20s : %4d (%4d)',files(i).name, ttl, total));
end;
disp('-----------------------------------------------------');
disp(sprintf('%0.0f files with %0.0f semicolons',totalI, total));