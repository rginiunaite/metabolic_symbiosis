clear all
% Implicit ODE solver

% fixed oxygen
oxygen1 = 0.5;
oxygen2 = 0.5;

IC(1) = 0;        % Initial extracellular lactate, 1st compartment
IC(2) = 0;        % Initial extracellular lactate, 2nd compartment
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;        % Initial population 2nd compartment
%IC(5) = 0.7;      % Initial oxygen 1st compartment    
%IC(6) = 0.3;      % Initial oxygen 2nd compartment  
IC(5) = 0.1;      % Initial MCT4 1st comparmtent
IC(6) = 0.1;      % Initial MCT4 2nd compartment
%IC(9)=0;
%IC(10)=0;

InitialPop = IC(3) ; % have variable for initial population in the first compartment

t0 = 0;
yp0 = [0 0 0 0 0 0 ]; % guess for initial values of derivatives
options=odeset('RelTol',1e-6);
[y0,yp0] = decic(@dynamics,t0,IC,[1 1 1 1 1 1 ],yp0,[0 0 0 0 0 0],options,InitialPop,oxygen1,oxygen2);
% first 4 components fixed, this is maximum, so that one could find consistent initial
% conditions for both derivative

T = 1e5;           % Sets the end of time interval
tspan = [0 T];
y0 = IC;           % Sets initial condition
options=odeset('RelTol',1e-4);

[t,y]=ode15i(@dynamics,tspan,y0,yp0,options,InitialPop,oxygen1,oxygen2);



% total population
figure

% Population
subplot(2,2,[1 2]);
plot(t,y(:,3)); % first population
hold on
plot(t,y(:,4)); % second population
plot(t,y(:,3) + y(:,4)); % total population
legend('1st component','2nd component','Total')
xlabel('time')
ylabel('Population')
title(['Populations with fixed oxygen c_1 = ' num2str(oxygen1) ' and c_2 = ' num2str(oxygen2) ' '])



% comparison of lactates
subplot(2,2,3);
plot(t,y(:,1));
hold on
plot(t,y(:,2));
legend('1st component','2nd component')
xlabel('time')
ylabel('Lactate')
title(['Lactate with fixed oxygen c_1 = ' num2str(oxygen1) ' and c_2 = ' num2str(oxygen2) ' '])


% comparison of MCT4
subplot(2,2,4);
plot(t,y(:,5));
hold on
plot(t,y(:,6));
legend('1st component','2nd component')
xlabel('time')
ylabel('MCT4')
title(['MCT4 with fixed oxygen c_1 = ' num2str(oxygen1) ' and c_2 = ' num2str(oxygen2) ' '])

function deriv = dynamics(t,y,yp,InitialPop,oxygen1, oxygen2)


    deriv=zeros(6,1); % create an empty column vector

    k1=1; % this is for intracellulr lactate
    k2=1; % this is for intracellulr lactate
    k5=0.005*InitialPop; % lactate degradation 
    k51=0.005*InitialPop;
    k52=0.5*0.005*InitialPop;

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;

    %omega12=0;
    %omega21=0;

    k4=0.001; % degradation of MCT4
    
    H0=1; % parameter that influences dynamics of MCT4, how fast it increases 
    %as a function of h_i

    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4

    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    
    if t<50000
        S2=1*InitialPop;    % in the hypoxic compartment source is the same
        %up to a certain time
    else
        S2=0;
        %  S2=1*InitialPop;
    end

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha
    
    if t<75000
        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
        h1=exp(beta1*(1-oxygen1));
        h2=exp(beta1*(1-oxygen2));
    else  
        k3=0.001;  
        %k3=0;
        h1=exp(beta1*(1-oxygen1));
        h2=exp(beta1*(1-oxygen2));
    end

    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   
% lactate 1

    deriv(1) = yp(1)-(-omega12*y(1) + omega21 * y(2) + y(3)* y(5) / (KMAX4 + y(1)) - ...
        y(3) * y (1)  / (KMAX1 + y(1)) - k5 * y(1)); % assume m1 = 1
    

% lactate 2      

    deriv(2) = yp(2)-(-omega21*y(2) + omega12 * y(1) + y(4)* y(6) / (KMAX4 + y(2)) - ...
        y(4) * y (2) / (KMAX1 + y(2)) - k5 * y(2));
 
% population 1

    deriv(3) = yp(3) - ((birth(oxygen1) - death(y(1),oxygen1)) * y (3));

% population 2

    deriv(4) = yp(4) - ((birth(oxygen2) - death(y(2),oxygen2)) * y (4));
    
% % oxygen 1
% 
%     deriv(5) = yp(5) - (S1 - diff12 * y(5) + diff21 * y(6) - k6 * y(3) * y(5));
%     
% % oxygen 2
%     
%     deriv(6) = yp(6) - (S2 - diff21 * y(6) + diff12 * y(5) - k6 * y(4) * y(6));
    
% MCT4 1

    deriv(5) = yp(5) - (k3 * h1 /(H0 + h1) - k4 * y (5));

% MCT4 2

    deriv(6) = yp(6) - (k3 * h2 /(H0 + h2) - k4 * y (6));

    
end


% birth of new cells dependent on oxygen
function value = birth (oxygen)

    % parameter values from de la Cruz et al. JTB (2016)

    exponent = -0.2;
    a0 = 8.25e3;
    oxycr = 0.02;

    value = a0 * ((oxygen/oxycr)-1)^exponent;
    value = 1/value;
end
%%%

% death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.1;

    value=nu0*lactate/(L0*oxygen+lactate);
end
%%%