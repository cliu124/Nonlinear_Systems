function G=loadG(fn)
%G=graph; 
ed=readtable(fn); ed(:,[1 2])
da=[ed.EndNodes_1(:),ed.EndNodes_2(:)]

ed2=table(da,'EndNodes'); ed2, pause 
G=graph(ed); 
return 
G.Edges=add
G.Nodes(:)=1:maxn; 