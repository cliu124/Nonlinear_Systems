function invtbl = inv_table(table_)
[nds,I] = sort(table_(:));
edgs = repmat(1:size(table_,1),1,size(table_,2))'; edgs = edgs(I);
invtbl = rpval2M(nds,edgs);