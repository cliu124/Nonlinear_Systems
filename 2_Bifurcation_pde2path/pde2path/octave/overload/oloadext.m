function p=oloadext(p); 
% extension of loadp for octave; p.pdeo not saved, hence needs to be regenerated
% in a dummy way; actual grid is saved in p.gr, and then restored. 
try
switch p.ndim;
  case 1; pde=stanpdeo1D(1,0.25); 
  case 2; pde=stanpdeo2D(1,1,1); 
  case 3; pde=stanpdeo3D(1,1,1,1); 
end
p.pdeo=pde; 
p.pdeo.grid.p=p.gr.po; p.pdeo.grid.t=p.gr.t; p.pdeo.grid.e=p.gr.e; 
p.gr=[]; 
end
