clear all;
close all;
clc;

%---------------
%using simulations with randon initial conditions. 
t_list=linspace(0,1000,10000); %simulation time, as 10 times the period
ran_num=100; %number of random initial conditions
ran_ini=randn(2,ran_num); %generate random initial conditions
for ran_index=1:ran_num %normalize the initial condition with unit norm
    ran_ini(:,ran_index)=ran_ini(:,ran_index)/norm(ran_ini(:,ran_index));
end
A=[-1/100,1;
    0,-1/10];

for t_ind=1:length(t_list)
    t=t_list(t_ind); %update time. 
    Phi=expm(A*t); %State transition matrix. 
    for ran_ind=1:ran_num
        x0=ran_ini(:,ran_ind); %select initial condition from random number liet. 
        x=Phi*x0; %evolve state variable x using state transition matrix Phi. 
        ran_run(ran_ind,t_ind)=norm(x)^2/norm(x0)^2; %get the energy growth from simulation results. 
    end
end

%-----------------
%compute upper bound of transient growth using Lyapunov method
%we solve: 
% min gamma,
% such that I<=P, for any t
% and P*A+A'*P+ dP/dt<=2*lambda*P, for any t. 

lambda=sdpvar(1,1); %upper bound of growth rate to be optimized
I=eye(2,2); %identity matrix
P=sdpvar(2,2);
constraint=[P>=I,P*A+A'*P<=2*lambda*P];

objective=lambda; %set objective of optimization as lambda. In default, it will minimize this objective function. If you want to maximize it, then set objective=-gamma
sdp_option=sdpsettings('solver','sedumi'); %This solver can be modified as mosek (https://www.mosek.com/), and using edu email, you can get a free license. 
results=bisection(constraint,objective,sdp_option); %call YALMIP function optimize to solve this SDP problem. This is a bi-linear problem, so it has to use bisection command instead of optimize!!!
%-------------------

%-------------------
%visulization of results. 
%the solution x(t) should be bounded as |x(t)|<= C*exp(\lambda*t);
for ran_ind=1:ran_num
    plot(t_list,log(sqrt(ran_run(ran_ind,:))),'--'); hold on; 
end
plot(t_list,(value(lambda)*t_list),'b--','Linewidth',3); hold on;
disp(['Growth rate from Lyapunov method: ',num2str(value(lambda))]);
disp(['Growth rate from eigenvalue: ',num2str(max(real(eig(A))))]);