clear all
% This is a verification that the population of two non-interacting species
% can be increased when they are different

a = 0.1; % parameter that makes the dynamics of two populations different 
%oxygen:
c1 = 0.5; 
c2 = 0.5;
% % initial populations
n0 = [10; 10];
tspan = [0 5];

% solve ode with fixed parameters
[t,y]=ode45(@(t,y) dn1dt(t,y,c1,c2,a),tspan,n0);

%calculate total population
total_population = y(:,1) + y(:,2);

figure
plot (t,total_population)
xlabel('time')
ylabel('n_1 + n_2')

c1 = 0.9; 
c2 = 0.1;
% solve ode with fixed parameters
[t,y]=ode45(@(t,y) dn1dt(t,y,c1,c2,a),tspan,n0);
%calculate total population
total_population_different = y(:,1) + y(:,2);

hold on
plot (t,total_population_different)
legend('c_1 = 0.5, c_2 =0.5','c_1 = 0.9, c_2 =0.1')
title('Total population for \alpha = 0.1')

% dynamics of the two populations
function value = dn1dt(t,n,c1,c2,a)

r1=1;
r2=1;
value = [ r1*n(1)*(c1-a*n(1)); r2*n(2)*(c2 - n(2))];
    
end
%%%
