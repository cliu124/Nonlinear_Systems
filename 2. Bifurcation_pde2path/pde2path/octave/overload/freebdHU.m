function es=freebdHU(x,t) % emulates delaunayTriangulation.freeBoundary; 
  % HU, 2020 
 if size(t,2)==3; % 2D, find edges which appear only once 
 e1=[t(:,[3 1]); t(:,[2,3]); t(:,[1 2])]; % all edges 
 e2=[t(:,[1 3]); t(:,[3,2]); t(:,[2 1])]; % inverse edges 
 si=~ismember(e1,e2,'rows'); es=e1(si,:); % find and extract singles 
else fprintf('3D\n'); % 3D, more options how faces of tetras can appear 
 e1=[t(:,[1,2,3]); t(:,[2,3,4]); t(:,[3,4,1]); t(:,[1,2,4])]; % all triangles in tetra-list 
 e2=[t(:,[1,3,2]); t(:,[2,4,3]); t(:,[3,1,4]); t(:,[1,4,2])]; % alternative orientations 
 e3=[t(:,[2,1,3]); t(:,[3,2,4]); t(:,[4,3,1]); t(:,[2,1,4])]; % 
 e4=[t(:,[2,3,1]); t(:,[3,4,2]); t(:,[4,1,3]); t(:,[2,4,1])]; % 
 e5=[t(:,[3,1,2]); t(:,[4,2,3]); t(:,[1,3,4]); t(:,[4,1,2])]; % 
 e6=[t(:,[3,2,1]); t(:,[4,3,2]); t(:,[1,4,3]); t(:,[4,2,1])]; % 
 ea=[e2;e3;e4;e5;e6];  
 si=~ismember(e1,ea,'rows'); es=e1(si,:); % remove doubles wrt e2,..,e6
 esl=size(es,1); es2=[]; essli=1:esl; % collect those that are also single in es 
 for i=1:esl; 
  cli=setdiff(essli,i); % the comparison list  
  if ~ismember(es(i,:),es(cli,:),'rows') es2=[es2; es(i,:)]; end 
 end 
 es=es2; 
 end 