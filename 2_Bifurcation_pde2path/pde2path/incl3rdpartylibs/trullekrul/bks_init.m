function bks = bks_init(tri)
nds = max(tri(:));
bks.rmnd = true(nds,1);
bks.mvnd = true(nds,1);