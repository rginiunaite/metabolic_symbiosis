clear all
% Implicit ODE solver
% Including intracellular lactate to reflect the effect of MCT correctly
% fixed MCT4 and MCT1 and they are turned on only when S2 = 0 
% with glucose dynamics



% note that I have to put here values of MCT which are inside the function
     alpha = 0.5;
     MCT42 = 1;
     MCT11 = MCT42;
     MCT12 = MCT42*alpha;
     MCT41 = MCT42*alpha;

%      m = 2;
%      beta = 0.5;
%      %eps = beta*m; % when relative to the increase
%      eps = 0;
%      MCT42 = m+eps;
%      MCT11 = m;
%      MCT12 = m;
%      MCT41 = m-eps;
          
% MCT11 = 1;
% MCT12 = 1;
% MCT41 = 0.6;
% MCT42 = 1.2;




IC(1) = 0;        % Initial extracellular lactate, 1st compartment
IC(2) = 0;        % Initial extracellular lactate, 2nd compartment
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;     % Initial population 2nd compartment
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
%IC(5) = 0.1;     % Initial MCT4 1st comparmtent
%IC(6) = 0.1;     % Initial MCT4 2nd compartment
IC(7) = 0;          % Intracellular lactate, 1st
IC(8) = 0;          % Intracellular lactate, 2nd
IC(9) = 0.1;      % Glucose, 1st
IC(10) = 0.1;     % Glucose, 2nd

InitialPop = IC(3) ; % have variable for initial population in the first compartment

t0 = 0;
yp0 = [0 0 0 0 0 0 0 0 0 0]; % guess for initial values of derivatives
options=odeset('RelTol',1e-6);
[y0,yp0] = decic(@dynamics,t0,IC,[1 1 1 1 1 1 1 1 1 1],yp0,[0 0 0 0 0 0 0 0 0 0],options,InitialPop,MCT41,MCT42,MCT11,MCT12);
% 1 corresponds to fixed components, 0 to variable, they are determined
% using this function

T = 2e6;           % Sets the end of time interval
tspan = [0 T];
y0 = IC;           % Sets initial condition
options=odeset('RelTol',1e-4);

[t,y]=ode15i(@dynamics,tspan,y0,yp0,options,InitialPop,MCT41,MCT42,MCT11,MCT12);

len = length(y(:,3)); % number of simulations done

final_normoxic = y(len,3);
final_hypoxic = y(len,4);
final_total = y(len,3) + y(len,4);

% total population
figure

% Population
subplot(3,3,[1 3]);
plot(t,y(:,3),'LineWidth', 2); % first population
hold on
plot(t,y(:,4),'LineWidth', 2); % second population
plot(t,y(:,3) + y(:,4),'LineWidth', 2); % total population
h_legend = legend('1st compartment','2nd compartment','Total');
set(h_legend,'FontSize',14)
%yaxis(0:2200)
xlabel('Time','FontSize',14)
ylabel('Population','FontSize',14)
%title(['Populations, final total ' num2str(final_total) ', final 1st component ' num2str(final_normoxic) ', final 2nd component ' num2str(final_hypoxic) ' '],'FontSize',14)
title({['Populations, final total ' num2str(final_total) ', final 1st compartment ' num2str(final_normoxic) ', final 2nd compartment ' num2str(final_hypoxic) ' '], ...
    [' MCT1_1 = ' num2str(MCT11) ', MCT1_2 = ' num2str(MCT12) ', MCT4_1 = ' num2str(MCT41), ' MCT4_2 = ' num2str(MCT42) ' ']},'FontSize',14)


% comparison of extracellular lactate
subplot(3,3,4);
plot(t,y(:,1),'LineWidth', 2);
hold on
plot(t,y(:,2),'LineWidth', 2);
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Lactate','FontSize',14)
title('Extracellular lactate ','FontSize',14)


% comparison of intracellular lactate
subplot(3,3,5);
plot(t,y(:,7),'LineWidth', 2);
hold on
plot(t,y(:,8),'LineWidth', 2);
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Lactate','FontSize',14)
title('Intracellular lactate ','FontSize',14)


% comparison of Oxygen
subplot(3,3,6);
plot(t,y(:,5),'LineWidth', 2);set(gca,'FontSize',14)
hold on
plot(t,y(:,6),'LineWidth', 2);set(gca,'FontSize',14)
%plot(t,y(:,5) + y(:,7)); % total population
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Oxygen','FontSize',14)
title('Oxygen','FontSize',14)


% comparison of Glucose
subplot(3,3,[7 9]);
plot(t,y(:,9),'LineWidth', 2);set(gca,'FontSize',14)
hold on
plot(t,y(:,10),'LineWidth', 2);set(gca,'FontSize',14)
%plot(t,y(:,5) + y(:,7)); % total population
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Glucose','FontSize',14)
title('Glucose','FontSize',14)


% figure
% plot([0 T],[MCT11 MCT11],[0 T],[MCT12 MCT12],[0 T],[MCT41 MCT41],[0 T],[MCT42 MCT42],'LineWidth',2)
% hlegend = legend('MCT11','MCT12','MCT41','MCT42')
% title([' MCT11 = ' num2str(MCT11) ', MCT12 ' num2str(MCT12) ', MCT41 ' num2str(MCT41), ' MCT42 ' num2str(MCT42) ' '],'FontSize',14)



% % comparison of MCT4
% subplot(2,2,4);
% plot(t,y(:,5),'LineWidth', 2);
% hold on
% plot(t,y(:,6),'LineWidth', 2);
% legend('1st component','2nd component')
% xlabel('time')
% ylabel('MCT4')
% title(['MCT4 with fixed oxygen c_1 = ' num2str(oxygen1) ' and c_2 = ' num2str(oxygen2) ' '])

function deriv = dynamics(t,y,yp,InitialPop,MCT41,MCT42,MCT11,MCT12)


    deriv=zeros(10,1); % create an empty column vector

    k1=1; % this is for intracellulr lactate
    k2=1; % this is for intracellulr lactate
    k5=0.005*InitialPop; % lactate degradation 
    k51=0.005;
    k52=0.5*0.005*InitialPop;

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1 = 0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4 = 0.1; % parameter for lactate dynamics, dependence on MCT4
   
    VM4 = 1;
    VM1 = 1;
    
    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component


    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = 0.5*BiasOxygenTransport*InitialPop;
   
% extracellular lactate 1

    deriv(1) = yp(1)-(-omega12*y(1) + omega21 * y(2) + k1*y(3)* MCT41*y(7)*VM4 / (KMAX4 + y(7)) - ...
        VM1*y(3) * y (1) * MCT11  / (KMAX1 + y(1)) - k5 * y(1)); % assume m1 = 1
    

% extracellular lactate 2      

    deriv(2) = yp(2)-(-omega21*y(2) + omega12 * y(1) + k1*y(4)* MCT42*y(8)*VM4 / (KMAX4 + y(8)) - ...
        VM1*y(4) * y (2) * MCT12 / (KMAX1 + y(2)) - k5 * y(2));
 
% population 1

    deriv(3) = yp(3) - ((birth1(y(5),y(9)) + birth2(y(9)) + birth3(y(5),y(7)) - death(y(7),y(5))) * y (3));

% population 2

    deriv(4) = yp(4) - ((birth1(y(6),y(10)) + birth2(y(10)) + birth3(y(6),y(8)) - death(y(8),y(6))) * y (4));
    
 % oxygen 1

    deriv(5) = yp(5) - (S1 - diff12 * y(5) + diff21 * y(6) - k6 * y(3) * y(5));
     
 % oxygen 2
     
    deriv(6) = yp(6) - (S2 - diff21 * y(6) + diff12 * y(5) - k6 * y(4) * y(6));
    
% intracellular lactate 

    % parameters for ODIL
    S01 = 50;
    S02 = 50;

    c0=0.05;
    
    % Glucose dependent source of intracellular lactate (GDIL)
    
    Sl1 = (S01*y(9))/(1+y(5)/c0); 

    Sl2 = (S02*y(10))/(1+y(6)/c0) ;
    

    deriv(7) = yp(7) - (Sl1 - VM4 * MCT41 * y(3)*y(7)/(KMAX4 + y(7)) + VM1*MCT11 * y(3)*y(1)/(KMAX1 +y(1)) - k51 *y(7)*y(3)); % 1st compartment
    
    deriv(8) = yp(8) - (Sl2 - VM4 * MCT42 * y(4)*y(8)/(KMAX4 + y(8)) + VM1*MCT12 * y(4)*y(2)/(KMAX1 +y(2)) - k51 * y(8)*y(4)); % 2nd compartment

    % Glucose
    Sg1 = 1000;
    Sg2 = 0;
    GL1 = 0;
    GL2 = 1000;
    k4 = 1;
    
    deriv(9) = yp(9) - ( Sg1 - GL1 * y(9) + GL2 * y(10) - k4 * y(9) * y(3));
    
    deriv(10) = yp(10) - ( Sg2 + GL1 * y(9) - GL2 * y(10) - k4 * y(10) * y(4));
     

end


% birth of new cells dependent on oxygen
function value = birth1 (oxygen,glucose)

    % parameter values from de la Cruz et al. JTB (2016)
    
    exponent = -0.2;
    a0 = 8.25e3;
    oxycr = 0.02^2;

    value = a0 * ((oxygen*glucose/oxycr)-1)^exponent;
    value = 1/value;
end
%%%

% birth of new cells dependent on oxygen
function value = birth2 (glucose)

    % parameter values from de la Cruz et al. JTB (2016)
    
    exponent = -0.2;
    a0 = 8.25e3*15;
    oxycr = 0.02;

    value = a0 * ((glucose/oxycr)-1)^exponent;
    value = 1/value;
end
%%%

% birth of new cells dependent on oxygen
function value = birth3 (oxygen,lactate)

    % parameter values from de la Cruz et al. JTB (2016)
    
    exponent = -0.2;
    a0 = 8.25e3*2;
    oxycr = 0.5*0.02^2;

    value = a0 * ((oxygen*lactate/oxycr)-1)^exponent;
    value = 1/value;
end
%%%


% death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.2; % previously was 0.1, to reduce the detrimental impact of lactate, do not knwo how it changed to 0.2 in his, or 0.1 in mine, I just wrote it wrong

    value=nu0*lactate/(L0*oxygen+lactate);
    
    value = 1e-3;
end
%%%