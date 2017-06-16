%clear all
% Implicit ODE solver
% fixed MCT4 
% only one variable for lactate deescribed in a very simple way





IC(1) = 0;        % Initial extracellular lactate, 1st compartment
IC(2) = 0;        % Initial extracellular lactate, 2nd compartment
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;        % Initial population 2nd compartment
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
%IC(5) = 0.1;      % Initial MCT4 1st comparmtent
%IC(6) = 0.1;      % Initial MCT4 2nd compartment
%IC(9)=0;
%IC(10)=0;

InitialPop = IC(3) ; % have variable for initial population in the first compartment



t0 = 0;
yp0 = [0 0 0 0 0 0]; % guess for initial values of derivatives
options=odeset('RelTol',1e-6);
[y0,yp0] = decic(@dynamics,t0,IC,[1 1 1 1 1 1],yp0,[0 0 0 0 0 0],options,InitialPop);
% 1 corresponds to fixed components, 0 to variable, they are determined
% using this function

T = 1e5;           % Sets the end of time interval
tspan = [0 T];
y0 = IC;           % Sets initial condition
options=odeset('RelTol',1e-4);

[t,y]=ode15i(@dynamics,tspan,y0,yp0,options,InitialPop);

len = length(y(:,3)); % number of simulations done


final_normoxic = y(len,3);
final_hypoxic = y(len,4);
final_total = y(len,3) + y(len,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Introducing radiation

% Radiation parameters

d = 4;
alpha = 0.35;
beta = 0.035;
m = 2.7;
K = 0.25;
OER_norm = 1;%(m*K + y(len,5))/(K + y(len,5)); % oxygen enhancement ratio
OER_hyp = 3;%(m*K + y(len,6))/(K + y(len,6));


%%%%%%
%Radiation%
%

final_normoxic_after_radiation = final_normoxic * exp(-alpha/OER_norm * d - beta/OER_norm^2 *d^2);
final_hypoxic_after_radiation = final_hypoxic * exp(-alpha/OER_hyp * d - beta/OER_hyp^2 *d^2);
final_total_after_radiation = final_normoxic_after_radiation + final_hypoxic_after_radiation;


figure
x = [final_normoxic final_normoxic_after_radiation; final_hypoxic final_hypoxic_after_radiation; final_total final_total_after_radiation];
%store = x;
xtl = {'Final normoxic' 'Final hypoxic' 'Final total'};
hb = bar(x);
set(gca, 'XtickLabel', xtl,'fontweight','bold','Fontsize',16)
 set(hb(1),'FaceColor','b')
 set(hb(2),'FaceColor','k')

ylabel('Population size','fontweight','bold','FontSize',16)
h_legend = legend('Before radiation','After radiation');
set(h_legend,'FontSize',16);



%%%% Figure, combination of two
figure 
xbar = [final_normoxic final_normoxic_after_radiation store(1,1) store(1,2);...
    final_hypoxic final_hypoxic_after_radiation store(2,1) store(2,2);...
    final_total final_total_after_radiation store(3,1) store(3,2)];
hbar = bar(xbar);
xtl = {'Final normoxic' 'Final hypoxic' 'Final total'};
set(gca, 'XtickLabel', xtl,'fontweight','bold','Fontsize',16)
 set(hbar(1),'FaceColor','b')
 set(hbar(2),'FaceColor','k')
 set(hbar(3),'FaceColor','y')
 set(hbar(4),'FaceColor',[0,0,0]+alpha)
 
 ylabel('Population size','fontweight','bold','FontSize',16)
h_legend = legend('Antiangiogenesis, before radiation','Antiangiogenesis, after radiation',...
    'Antiangiogenesis + MCT4 downregulation, before radiation','Antiangiogenesis + MCT4 downregulation, after radiation');
set(h_legend,'FontSize',16);
 set(gca,'linewidth',2)

% figure 
% plot(mct,log(final_normoxic),mct,log(final_hypoxic),mct,log(final_total),'LineWidth', 2)
% h_legend = legend('1st component','2nd component','Total')
% set(h_legend,'FontSize',14)
% xlabel('MCT','FontSize',14)
% ylabel('log(Population)','FontSize',14)
% title(['MCT influence on population sizes, final total ' num2str(final_total(i)) ', final normoxic ' num2str(final_normoxic(i)) ', final hypoxic ' num2str(final_hypoxic(i)) ' '],'FontSize',14)
% set(gca,'XTick',[1 2 3 4 5 6] ); %This are going to be the only values affected.
% set(gca,'XTickLabel',[10^-5 10^-4 10^-3 10^-2 10^-1 1] ); %This is what it's going to appear in those places.
% grid on
% 




% total population
figure

% Population
subplot(3,1,1);
plot(t,y(:,3),'LineWidth', 2); % first population
hold on
plot(t,y(:,4),'LineWidth', 2); % second population
plot(t,y(:,3) + y(:,4),'LineWidth', 2); % total population
h_legend = legend('1st compartment','2nd compartment','Total');
set(h_legend,'FontSize',14);
xlabel('Time','fontweight','bold','FontSize',14)
ylabel('Population','fontweight','bold','FontSize',14)
title(['Populations, final total ' num2str(final_total) ', final 1st component ' num2str(final_normoxic) ', final 2nd component ' num2str(final_hypoxic) ' '],'FontSize',14)
set(gca,'linewidth',2)


% comparison of lactates
subplot(3,1,2);
plot(t,y(:,1),'LineWidth', 2);
hold on
plot(t,y(:,2),'LineWidth', 2);
% h_legend = legend('1st component','2nd component')
% set(h_legend,'FontSize',14)
xlabel('Time','fontweight','bold','FontSize',14)
ylabel('Lactate','fontweight','bold','FontSize',14)
title('Lactate ','FontSize',14)
set(gca,'linewidth',2)

% comparison of Oxygen
subplot(3,1,3);
plot(t,y(:,5),'LineWidth', 2);%set(gca,'FontSize',14)
hold on
plot(t,y(:,6),'LineWidth', 2);%set(gca,'FontSize',14)
%plot(t,y(:,5) + y(:,7)); % total population
% h_legend = legend('1st component','2nd component');
% set(h_legend,'FontSize',14)
xlabel('Time','fontweight','bold','FontSize',14)
ylabel('Oxygen','fontweight','bold','FontSize',14)
title('Oxygen','FontSize',14)
set(gca,'linewidth',2)

% % comparison of MCT4
% subplot(2,2,4);
% plot(t,y(:,5),'LineWidth', 2);
% hold on
% plot(t,y(:,6),'LineWidth', 2);
% legend('1st component','2nd component')
% xlabel('time')
% ylabel('MCT4')
% title(['MCT4 with fixed oxygen c_1 = ' num2str(oxygen1) ' and c_2 = ' num2str(oxygen2) ' '])

function deriv = dynamics(t,y,yp,InitialPop)


    deriv=zeros(6,1); % create an empty column vector

    k1=1; % this is for intracellulr lactate
    k2=1; % this is for intracellulr lactate
    k5=0.005*InitialPop; % lactate degradation 
    k51=0.005*InitialPop;
    k52=0.5*0.005*InitialPop;

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=100;%0;0.1*InitialPop; % transporters
    omega21=500;%BiasLactateTransport*InitialPop;

%     omega12=0;
%     omega21=0;


    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component
%     if t<50000
%         S2=1*InitialPop;    % in the hypoxic compartment source is the same
%         %up to a certain time
%     else
%         S2=0;
%         S2=1*InitialPop;
%     end
    

    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 100;%0.05*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = 100;%BiasOxygenTransport*0.05*InitialPop;
   
% lactate 1


    % add source
    S01=10;
    S02=10;

    c0=0.05;
    
    % Oxygen dependent source of intracellular lactate (ODIL)
    
    Sl1=S01/(1+(y(5)/c0)); 

    Sl2=S02/(1+(y(6)/c0));


    deriv(1) = yp(1)-(Sl1 -omega12*y(1) + omega21 * y(2) - k5 * y(1)); % assume m1 = 1
    

% lactate 2      

    deriv(2) = yp(2)-(Sl2 -omega21*y(2) + omega12 * y(1) - k5 * y(2));
 
% population 1

    deriv(3) = yp(3) - ((birth(y(5)) - death(y(1),y(5))) * y (3));

% population 2

    deriv(4) = yp(4) - ((birth(y(6)) - death(y(2),y(6))) * y (4));
    
 % oxygen 1

    deriv(5) = yp(5) - (S1 - diff12 * y(5) + diff21 * y(6) - k6 * y(3) * y(5));
     
 % oxygen 2
     
     deriv(6) = yp(6) - (S2 - diff21 * y(6) + diff12 * y(5) - k6 * y(4) * y(6));
    
% % MCT4 1
% 
%     deriv(5) = yp(5) - (k3 * h1 /(H0 + h1) - k4 * y (5));
% 
% % MCT4 2
% 
%     deriv(6) = yp(6) - (k3 * h2 /(H0 + h2) - k4 * y (6));

    
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
    L0 = 0.2;

    value=nu0*lactate/(L0*oxygen+lactate);
end
%%%