% Checking symbiosis for a very simplistic model


tspan = [0 5];

y0 =[5;5;10;10]; % initially lactate concentration is 5 in both components 
% population is 10

% solve ode with fixed parameters
[t,y]=ode45( @deriv,tspan,y0);

figure
plot(t,y(:,1),t,y(:,3)+y(4))

% lactate dynamics

function value = deriv(t,y)

% Oxygen concentration

c1 = 0.5;
c2 = 0.5;


% Lactate transporters MCT4 anad MCT1, assume that they are fixed
mct4 = [1;10]; % different transporters are activated in different components
%mct4 = [0;0]; % transporters are off

% other way arounf for MCT1 transporter
mct1 = [10,1]; % different transporters are activated in different 
%components, mutually exclusive from MCT4
%mct1 = [0;0]; % transporters are off


% parameters for lactate dynamics
w = [1;1]; % diffusion
VM4 = 1;
KM4 = 1;
VM1 = 1;
KM1 = 1;
k5 = 1;

% parameters for cell population dynamics

a0 = 5;
oc = 0.8;
pow = -1;
v0 = 1;


% lactate 1
    value(1) = -w(1)*y(1) + w(2) * y(2) + VM4/(KM4 + y(1))*mct4(1)*y(3) - ...
        VM1/(KM1 +y(1)) * mct1(1)*y(3)*y(1) - k5 *y(1);
        
% lactate 2      
    value(2) = -w(2)*y(2) + w(1) * y(1) + VM4/(KM4 + y(2))*mct4(2)*y(4) - ...
        VM1/(KM1 +y(2)) * mct1(2)*y(4)*y(2) - k5 *y(2);
    
% population 1

    value(3) = (a0 * (c1/oc - 1)^pow - v0 * y(1)/(c1 + y(1))) * y(3); % assume L(c) = c;

% population 2

    value(4) = (a0 * (c2/oc - 1)^pow - v0 * y(2)/(c1 + y(2))) * y(4); % assume L(c) = c;  
    
    value = value'; % since column vector needs to be returned
end


