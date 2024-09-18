This repository is example code for 
Nonlinear Systems ME 3295-001/ME 5895-001/ECE 6095-004
at University of Connecticut, taught by Dr. Chang Liu (https://changliulab.engineering.uconn.edu/).

This provide examples to use YALMIP (https://yalmip.github.io/) for searching Lyapunov function and pde2path (https://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/) to conduct bifurcation analysis.

These softwares have been already downloaded and unzipped. 

Running install.m will finish the installation of YALMIP, SeDuMi and pde2path.

For YALMIP, this repository already has SeDuMi (https://sedumi.ie.lehigh.edu/?page_id=58) as semi-definite programming solver.

The Mosek (https://www.mosek.com/) solver is required for A_growth_rate.m, C_region_of_attraction.m, D_sum_of_squares.m, E_growth_rate_time_varying.m. This is NOT contained here and it needs to install Mosek following the instructions in the link (https://docs.mosek.com/10.2/install/installation.html). 

Using edu email can get a free academic license of Mosek. 