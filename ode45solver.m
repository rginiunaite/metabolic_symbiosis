clear all
% ODE45 solver for the system


oxygen1 = 0.9;
oxygen2 = 0.1;

IC(1) = 0;        % Initial extracellular lactate, 1st compartment
IC(2) = 0;        % Initial extracellular lactate, 2nd compartment
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;        % Initial population 2nd compartment
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
IC(7) = 0.1;      % Initial MCT4 1st comparmtent
IC(8) = 0.1;      % Initial MCT4 2nd compartment
%IC(9)=0;
%IC(10)=0;

InitialPop = IC(3); % have variable for initial population in the first compartment


T = 100;           % Sets the end of time interval
tspan = [0 T];

[t,y]=ode45(@(t,y) dynamics(t,y,InitialPop),tspan,IC); 

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
% subplot(6,1,1);
% plot(notrt,notry(:,3)+notry(:,4));
% hold on
% plot(t,y(:,3)+y(:,4));
% legend('No lactate transporters','lactate transporters')
% xlabel('time')
% ylabel('n_1 + n_2')
% title(['Total population, initially ' num2str(y0(3)) ' in each of the components'])
% % 
% %lactate in total
% subplot(6,1,2);
% plot(notrt,notry(:,1)+notry(:,2));
% hold on
% plot(t,y(:,1)+y(:,2));
% legend('No lactate transporters','lactate transporters')
% xlabel('time')
% ylabel('l_1 + l_2')
% title('Total lactate')

% % no transporters, comparison of the two populations
% subplot(6,1,3);
% plot(notrt,notry(:,3));
% hold on
% plot(notrt,notry(:,4));
% legend('1st component','2nd component')
% xlabel('time')
% ylabel('Population')
% title('Populations no transporters')

% with transporters, comparison of the two populations
subplot(2,1,1);
plot(t,y(:,3));
hold on
plot(t,y(:,4));
legend('1st component','2nd component')
xlabel('time')
ylabel('Population')
title('Populations with transporters')

% %no transporters, comparison of lactates
% subplot(6,1,5);
% plot(notrt,notry(:,1));
% hold on
% plot(notrt,notry(:,2));
% legend('1st component','2nd component')
% xlabel('time')
% ylabel('Lactate')
% title('Lactate no transporters')

%with transporters, comparison of lactates
subplot(2,1,2);
plot(t,y(:,1));
hold on
plot(t,y(:,2));
legend('1st component','2nd component')
xlabel('time')
ylabel('Lactate')
title('Lactate with transporters')



%function deriv = dynamics(t,y,yp,InitialPop)
function deriv = dynamics(t,y,InitialPop)

    deriv=zeros(8,1); % create an empty column vector

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

    beta1=2.5; %parameter for the hypoxia dependent on HIF-1alpha
    
    if t<75000
        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
        h1=exp(beta1*(1-y(5)));
        h2=exp(beta1*(1-y(6)));
    else  
        %k3=0.001;  
        k3=0;
        h1=exp(beta1*(1-y(5)));
        h2=exp(beta1*(1-y(6)));
    end

    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   

    % lactate 1

    deriv(1) = (-omega12*y(1) + omega21 * y(2) + y(3)* y(7) / (KMAX4 + y(1)) - ...
        y(3) * y (1) / (KMAX1 + y(1)) - k5 * y(1));
        
% lactate 2      

    deriv(2) = (-omega21*y(2) + omega12 * y(1) + y(4)* y(8) / (KMAX4 + y(2)) - ...
        y(4) * y (2) / (KMAX1 + y(2)) - k5 * y(2));
 
% population 1

    deriv(3) =  ((birth(y(5)) - death(y(1),y(5))) * y (3));

% population 2

    deriv(4) = ((birth(y(6)) - death(y(2),y(6))) * y (4));
    
% oxygen 1

    deriv(5) =  (S1 - diff12 * y(5) + diff21 * y(6) - k6 * y(3) * y(5));
    
% oxygen 2
    
    deriv(6) = (S2 - diff21 * y(6) + diff12 * y(5) - k6 * y(4) * y(6));
    
% MCT4 1

    deriv(7) =  (k3 * h1 /(H0 + h1) - k4 * y (7));

% MCT4 2

    deriv(8) =  (k3 * h2 /(H0 + h2) - k4 * y (8));
    
    
    
end


% birth of new cells dependent on oxygen
function value = birth (oxygen)

    % parameter values from de la Cruz et al. JTB (2016)

    exponent = -0.2;
    a0 = 8.25e3;
    oxycr=0.02;

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
