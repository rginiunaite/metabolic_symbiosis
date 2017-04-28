clear all
% Checking symbiosis for a very simplistic model


tspan = [0 10];

y0 =[0.5;0.5;100;100]; % initially lactate concentration is 5 in both components 
% population is 10

% Lactate transporters MCT4 anad MCT1, assume that they are fixed
mct4 = [1;10]; % different transporters are activated in different components
nomct4 = [0;0]; % transporters are off

% other way arounf for MCT1 transporter
mct1 = [10,1]; % different transporters are activated in different 
%components, mutually exclusive from MCT4
nomct1 = [0;0]; % transporters are off


% solve ode with fixed parameters
[t,y]=ode45(@(t,y) deriv(t,y,nomct4,nomct1),tspan,y0);

notrt = t; %no transporters t
notry = y; % no transporters y


[t,y]=ode45(@(t,y) deriv(t,y,mct4,mct1),tspan,y0); %adding transporters

% hold on 
% subplot(2,1,1);
% plot(t,ytr(:,3)+ytr(:,4))
% 
% subplot(2,1,2);
% plot(t,ytr(:,3)+ytr(:,4))
% 
% xlabel('time')
% ylabel('n_1 + n_2')
% 
% legend('No lactate transporters','lactate transporters')
% %title('Total population for \alpha = 0.1')


% total population
figure
subplot(2,1,1);
plot(notrt,notry(:,3)+notry(:,4));
hold on
plot(t,y(:,3)+y(:,4));
legend('No lactate transporters','lactate transporters')
xlabel('time')
ylabel('n_1 + n_2')
title('Populations')


subplot(2,1,2);
plot(notrt,notry(:,1)+notry(:,2));
hold on
plot(t,y(:,1)+y(:,2));
legend('No lactate transporters','lactate transporters')
xlabel('time')
ylabel('l_1 + l_2')
title('Lactate')



% lactate dynamics

function value = deriv(t,y,mct4,mct1)

% Oxygen concentration

c1 = 0.5;
c2 = 0.5;



% parameters for lactate dynamics
w = [0.5;0.5]; % diffusion
VM4 = 1;
KM4 = 1;
VM1 = 1;
KM1 = 1;
k5 = 0.1;

% parameters for cell population dynamics

a0 = 1;
oc = 1;
pow = 1;
v0 = 3;


% lactate 1
    value(1) = -w(1)*y(1) + w(2) * y(2) - VM4/(KM4 + y(1))*mct4(1)*y(3)*y(1)*y(1) + ...
        VM1/(KM1 +y(1)) * mct1(1)*y(3)*y(2) - k5 *y(1);
        
% lactate 2      
    value(2) = -w(2)*y(2) + w(1) * y(1) - VM4/(KM4 + y(2))*mct4(2)*y(4)*y(2) + ...
        VM1/(KM1 +y(2)) * mct1(2)*y(4)*y(1) - k5 *y(2);
    
% other way of defining dynamics of lactate, no constants but should be
% added 
   
% % lactate 1
%     value(1) = -w(1)*y(1) + w(2) * y(2) + mct1(1) * y(2)/(y(2)+y(1)) - ...
%         mct4(1)*y(1) - k5 *y(1);
%         
% % lactate 2      
%     value(2) = -w(2)*y(2) + w(1) * y(1) + mct1(2) * y(1)/(y(1)+y(2)) - ...
%         mct4(2)*y(2) - k5 *y(2);
 
% population 1

    value(3) = (a0 * (c1/oc )^pow - v0 * y(1)/(c1 + y(1))) * y(3); % assume L(c) = c;

% population 2

    value(4) = (a0 * (c2/oc )^pow - v0 * y(2)/(c1 + y(2))) * y(4); % assume L(c) = c;  
    
    value = value'; % since column vector needs to be returned
end


