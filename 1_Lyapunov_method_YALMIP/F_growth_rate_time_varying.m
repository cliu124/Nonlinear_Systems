clear all;
close all;
clc;

%---------------
%initialize the time-varying A matrix, and compute the growth rate
%We consider Example 4.22 in Khalil's textbook. 

%using simulations with randon initial conditions. 
sample_num=500; %This is the sampling number over time for one period
period_num=10; %number of sampling period
t_list=linspace(0,period_num*2*pi,period_num*sample_num); %simulation time, as 10 times the period
ran_num=100; %number of random initial conditions
ran_ini=randn(2,ran_num); %generate random initial conditions
for ran_index=1:ran_num %normalize the initial condition with unit norm
    ran_ini(:,ran_index)=ran_ini(:,ran_index)/norm(ran_ini(:,ran_index));
end

dt=diff(t_list); dt=dt(1); %get dt, assuming uniform space in time. 

Phi=eye(2,2); %state transition matrix at t=0. 
for t_ind=1:length(t_list)
    t=t_list(t_ind); %update time. 
    A(:,:,t_ind)=[-1+1.5*cos(t)^2, 1-1.5*sin(t)*cos(t);
       -1-1.5*sin(t)*cos(t), -1+1.5*sin(t)^2]; %A matrix for a time-varying system of Example 4.22 in Khalil's textbook
    Phi=expm(A(:,:,t_ind)*dt)*Phi; %update state transition matrix. 
    for ran_ind=1:ran_num
        x0=ran_ini(:,ran_ind); %select initial condition from random number liet. 
        x=Phi*x0; %evolve state variable x using state transition matrix Phi. 
        ran_run(ran_ind,t_ind)=norm(x)^2/norm(x0)^2; %get the transient growth from simulation results. 
    end
end

%-----------------
%compute upper bound of transient growth using Lyapunov method
%we solve: 
% min gamma,
% such that I<=P(t), for any t
% and P(t)A(t)+A(t)'P(t)+ dP/dt<=2*lambda*P(t), for any t. 

lambda=sdpvar(1,1); %upper bound of growth rate to be optimized
I=eye(2,2); %identity matrix
constraint=[]; %constraint list. 
for t_ind=1:sample_num% we only need to go through one period as this A is periodic.
    P{t_ind}=sdpvar(2,2); %a time-varying Lyapunov function
    constraint=[constraint,P{t_ind}>=I]; %add constraint that I<=P(t), for any t
end
P{sample_num+1}=P{1}; %enforce that P is also periodic. This is only suitable for periodic A(t). 

for t_ind=1:sample_num
    dP=(P{t_ind+1}-P{t_ind})/(dt); %use finite difference to approximate  dP/dt
    constraint=[constraint,P{t_ind}*A(:,:,t_ind)+A(:,:,t_ind)'*P{t_ind}+dP<=2*lambda*P{t_ind}]; %add Lyapunov inequality constraint: P*A+A'*P+dP/dt<=2\lambda*P(t), for any t. 
end

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