clear all
% Implicit ODE solver
% Including intracellular lactate to reflect the effect of MCT correctly
% fixed MCT4 and MCT1 and they are turned on only when S2 = 0 
% with glucose dynamics, simplified version



% note that I have to put here values of MCT which are inside the function
     alpha = 0.2;
     MCT42 = 1;
     MCT11 = MCT42;
     MCT12 = MCT42*alpha;
     MCT41 = MCT42*alpha;



IC(1) = 0;          % lactate, 1st
IC(2) = 0;          % lactate, 2nd
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;     % Initial population 2nd compartment
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
%IC(5) = 0.1;     % Initial MCT4 1st comparmtent
%IC(6) = 0.1;     % Initial MCT4 2nd compartment

IC(7) = 0.1;      % Glucose, 1st
IC(8) = 0.1;     % Glucose, 2nd

InitialPop = IC(3) ; % have variable for initial population in the first compartment

t0 = 0;
yp0 = [0 0 0 0 0 0 0 0 ]; % guess for initial values of derivatives
options=odeset('RelTol',1e-6);
[y0,yp0] = decic(@dynamics,t0,IC,[1 1 1 1 1 1 1 1 ],yp0,[0 0 0 0 0 0 0 0 ],options,InitialPop,MCT41,MCT42,MCT11,MCT12);
% 1 corresponds to fixed components, 0 to variable, they are determined
% using this function

T = 2e5;           % Sets the end of time interval
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
subplot(2,3,[1 3]);
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
title({['Populations, final total ' num2str(real(final_total)) ', final 1st compartment ' num2str(real(final_normoxic)) ', final 2nd compartment ' num2str(real(final_hypoxic)) ' '], ...
    [' MCT1_1 = ' num2str(MCT11) ', MCT1_2 = ' num2str(MCT12) ', MCT4_1 = ' num2str(MCT41), ' MCT4_2 = ' num2str(MCT42) ' ']},'FontSize',14)


% comparison of extracellular lactate
subplot(2,3,4);
plot(t,y(:,1),'LineWidth', 2);
hold on
plot(t,y(:,2),'LineWidth', 2);
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Lactate','FontSize',14)
title('Lactate ','FontSize',14)


% comparison of intracellular lactate
subplot(2,3,5);
plot(t,y(:,7),'LineWidth', 2);
hold on
plot(t,y(:,8),'LineWidth', 2);
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Glucose','FontSize',14)
title('Glucose','FontSize',14)


% comparison of Oxygen
subplot(2,3,6);
plot(t,y(:,5),'LineWidth', 2);%set(gca,'FontSize',14)
hold on
plot(t,y(:,6),'LineWidth', 2);%set(gca,'FontSize',14)
%plot(t,y(:,5) + y(:,7)); % total population
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Oxygen','FontSize',14)
title('Oxygen','FontSize',14)


function deriv = dynamics(t,y,yp,InitialPop,MCT41,MCT42,MCT11,MCT12)


    deriv=zeros(8,1); % create an empty column vector


    k51 = 1; % lactate consumption per cell
    k6 = 1; % oxygen consumption per cell
    k4 = 1; % glucose consumption per cell

        
    S1 = 2; %0.005*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component

   
    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp, oxygen diffusion
    diff12 = 0.1;
    diff21 = 0.5*BiasOxygenTransport*InitialPop;
    diff21 = 0.1;
    %  lactate 

%     % parameters for ODIL
%     S01 = 10;
%     S02 = 10;
% 
%     c0=0.05;
%     
%     % Glucose dependent source of intracellular lactate (GDIL)
%     
%     Sl1 = (S01*y(7))/(1+y(5)/c0); 
% 
%     Sl2 = (S02*y(8))/(1+y(6)/c0) ;
    
    
    KO = 2e-4;
    gG = 1;
    gL = 1;
    gO = 1;
    
    % Glucose parameters
    
    Sg1 = 1;
    Sg2 = 0;
    GL1 = 0.1; % glucose diffusion from the first to the second compartment
    GL2 = 0; % glucose diffusion from the second to the first compartment
   
    KG = 1;
    gamL = 1;
    
    
    gamG = 1;
    BL = 4.63e-4;
    BG = 2.78e-3;
    BO = 2.32e-4; 
    k41 = 0.01; % lactate production as a ruslt of glycolysis
    
    % the rate of lactate consumption 
    
    kl1 = BL * y(5)/(y(5) + KO)* (gL * y(1))/(gG*y(7)+gL * y(1) + gamL);
    kl2 = BL *  y(6)/(y(6) + KO)* (gL * y(2))/(gG*y(8)+gL * y(2) + gamL);
    
    % the rate of glucose consumption
    
    kg1 = BG * y(7)/(y(7)+ gO * y(5)+KG);
    kg2 = BG * y(8)/(y(8)+ gO * y(6)+KG);
    
    % the rate of oxygen consumption
    
    ko1 = BO * y(5)/(y(5)+KO) * (gG*y(7)/ (gG*y(7)+ gL *y(1) + gamG));
    ko2 = BO * y(6)/(y(6)+KO) * (gG*y(8)/ (gG*y(8)+ gL *y(2) + gamG));
    

    deriv(1) = yp(1) - ( - MCT41 * y(1) + MCT11 * y(2) - k51 * kl1 * y(3) + 2 * y (3) *k41 * kg1);
    
    deriv(2) = yp(2) - ( - MCT42 *y(2) + MCT12 * y(1) - k51 * y(4) * kl2 + 2 * y(4) * k41 * kg2); % 2nd compartment
    
 
% population 1

    deriv(3) = yp(3) - ((birth1(y(5),y(7)) + birth2(y(7)) + birth3(y(5),y(1)) - death(y(1),y(5))) * y (3));

% population 2

    deriv(4) = yp(4) - ((birth1(y(6),y(8)) + birth2(y(8)) + birth3(y(6),y(2)) - death(y(2),y(6))) * y (4));
    
 % oxygen 1

    deriv(5) = yp(5) - (S1 - diff12 * y(5) + diff21 * y(6) - 6 * k6 * y(3) * ko1  - 3 * k51* y(3) * kl1);
     
 % oxygen 2
     
    deriv(6) = yp(6) - (S2 - diff21 * y(6) + diff12 * y(5) - 6 * k6 * y(4) * ko2 - 3* k51 * y(4) * kl2);
    


    % Glucose

    
    deriv(7) = yp(7) - ( Sg1 - GL1 * y(7) + GL2 * y(8) - k4 *  y(3) * kg1 - k6 * y(3) * ko1);
    
    deriv(8) = yp(8) - ( Sg2 + GL1 * y(7) - GL2 * y(8) - k4 *  y(4) * kg2 - k6 * y(4) * ko2);
     

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
    oxycr = 0.5*0.02^2; % since I need half less oxygen molecules for a reaction

    value = a0 * ((oxygen^1.5*lactate/oxycr)-1)^exponent; % maybe changed oxygen power to 1.5 to ensure that lactate is much more beneficial if there is a lot of oxygen
    value = 1/value;
end
%%%


% death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.2; 

    value=nu0*lactate/(L0*oxygen+lactate) + 10e-4;
    
 %   value = 1e-3;
end
%%%