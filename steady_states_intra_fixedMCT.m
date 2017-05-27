clear all
% this is to find steady states numerically, for dynamcis with
% intracellular lactate
% a little bit approximate value
% fixed MCT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
    % fixed oxygen (temporal)
    


    
    InitialPop =1000;

    k1=1; % this is for intracellulr lactate
    k2=1; % this is for intracellulr lactate
    k5=0.005*InitialPop; % lactate degradation 
    k51=0.005*InitialPop;
    k52=0.5*0.005*InitialPop;

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;

%     omega12=0;
%     omega21=0;

    k4=0.001; % degradation of MCT4
    
    H0=1; % parameter that influences dynamics of MCT4, how fast it increases 
    %as a function of h_i

  

    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component
    

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    

        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia



    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




IC(1) = 0.1;       % Initial extracellular lactate, 1st compartment
IC(2) = 0.1;       % Initial extracellular lactate, 2nd compartment
% choose carefully
% for KMAX need something as big as 6000 500
IC(3) = 5000;         % Initial population 1st compartment % similar to solution
IC(4) = 500;         % Initial population 2nd compartment % similar to solution
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
%IC(5) = 0.1;       % Initial MCT4 1st comparmtent
%(6) = 0.1;       % Initial MCT4 2nd compartment
IC(7)=0.1;           % Initial extracellular lactate, 1st compartment
IC(8)=0.1;           % Initial extracellular lactate, 2nd compartment


InitialPop = IC(3); % have variable for initial population in the first compartment

count = 0;
MCT11 = 1;
MCT12 = 1;
MCT41 = 0.6; % optimal for general case, consider what happens
MCT42 = 1.2;


VM4 = 1;
VM1 = 1;
KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
KMAX4=1; % parameter for lactate dynamics, dependence on MCT4

NumberOfIterations = 100;

% initialising final quantities to increase speed
final_normoxic = zeros(1,NumberOfIterations);
final_hypoxic = zeros(1,NumberOfIterations);
final_total = zeros(1,NumberOfIterations);
final_normoxic_EL = zeros(1,NumberOfIterations);
final_hypoxic_EL = zeros(1,NumberOfIterations);
final_normoxic_IL = zeros(1,NumberOfIterations);
final_hypoxic_IL = zeros(1,NumberOfIterations);
final_normoxic_O = zeros(1,NumberOfIterations);
final_hypoxic_O = zeros(1,NumberOfIterations);


for i = 0:NumberOfIterations %%% this is for VM1[0.5 1 10 100]
    
    count = count +1;
    
%   % varying KMAX4
%   KMAX4 = 0.1*i;
%     
%   % instead of zero take something very small
%     if i == 0
%        KMAX4 = 0.05;
%     end
    
   
      % varying KMAX1
  %   KMAX1 = 0.1*i; % this was not good at all
    
%   % instead of zero take something very small
%     if i == 0
%        KMAX1 = 0.05;
%    end


    
   MCT11 = 0.01 * i; % varying MCT1 as well in mutually exclusive regions
   MCT12 = 0.005 * i;
   
   %VM4 = 0.01*i;
   
   % VM1 = i;

    MCT41 = 0.005 * i;
    MCT42 = 0.01 * i;
    
   
    % for the size of MCT4 in one compartment relative to the other (i=0:5)
    %MCT42 = 0.1*(i+3); % 
    %MCT41 = 0.5;
    
    fun = @(y) system(y,MCT41, MCT42, MCT11, MCT12, VM4, KMAX4,VM1,KMAX1);
   
    options=optimset('disp','iter','LargeScale','off','TolFun',1e-7,'MaxIter',1e8,'MaxFunEvals',1e8);

   [y,fval,exitflag,output] = fsolve(fun,IC,options);

    len = length(y(:,3)); % number of simulations done

    final_normoxic(count) = y(len,3);
    final_hypoxic(count) = y(len,4);
    final_total(count) = y(len,3) + y(len,4);
    
    final_normoxic_EL(count) = y(len,1);
    final_hypoxic_EL(count) = y(len,2);
    %final_total_EL(count) = y(len,1) + y(len,2);
    
    final_normoxic_IL(count) = y(len,7);
    final_hypoxic_IL(count) = y(len,8);
    %final_total(count) = y(len,3) + y(len,4);
    
    final_normoxic_O(count) = y(len,5);
    final_hypoxic_O(count) = y(len,6);
    %final_total(count) = y(len,3) + y(len,4);

   
end
   

mct=0:NumberOfIterations; % defines how many iterations were made
%mct =0:3;

figure 
subplot(4,1,1)
plot(mct,(final_normoxic),mct,(final_hypoxic),mct,(final_total),'LineWidth', 2)
h_legend = legend('1st compartment','2nd compartment','Total');
set(h_legend,'FontSize',14)
xlabel('MCT4_2','FontSize',14)
%xlabel('V_{M_4}','FontSize',14)
%xlabel('K_{M_4}','FontSize',14)
%xlabel('V_{M_1}','FontSize',14)

ylabel('St. st. population','FontSize',14)
%title(['MCT4 influence on population sizes, MCT1 = ', num2str(MCT11)],'FontSize',14)
title('MCT4 and MCT1 influence on population sizes','FontSize',14)
%title('V_{M_4} influence on population sizes','FontSize',14)
%title('K_{M_4} influence on population sizes','FontSize',14)
%title('V_{M_1} influence on population sizes','FontSize',14)

%This are going to be the only values affected for the case when varying MCT4
set(gca,'XTick',[0:20:400] ); 
set(gca,'XTickLabel',[0:0.2:4] ); 
grid on
% % for varying VM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.01:NumberOfIterations*0.01] );
% grid on

% % for varying KM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.1:NumberOfIterations*0.1] );
% grid on

% % for varying VM1
% set(gca,'XTick',[0:4] );
% set(gca,'XTickLabel',[0.1 1 10 100] );
% grid on


subplot(4,1,2)
plot(mct,(final_normoxic_EL),mct,(final_hypoxic_EL),'LineWidth', 2)
% h_legend = legend('1st component','2nd component','Total');
% set(h_legend,'FontSize',14)
xlabel('MCT4_2','FontSize',14)
%xlabel('V_{M_4}','FontSize',14)
%xlabel('K_{M_4}','FontSize',14)
%xlabel('V_{M_1}','FontSize',14)

ylabel('St. st. extra lac','FontSize',14)
title('MCT4 and MCT1 influence on extracellular lactate','FontSize',14)
%title('V_{M_4} influence on extracellular lactate','FontSize',14)
%title('K_{M_4} influence on extracellular lactate','FontSize',14)
%title('V_{M_1} influence on extracellular lactate','FontSize',14)

%MCT4 axis
set(gca,'XTick',[0:20:400] ); %This are going to be the only values affected.
set(gca,'XTickLabel',[0:0.2:4] ); %This is what it's going to appear in those places.
grid on
% % for varying VM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.01:NumberOfIterations*0.01] );


% % for varying KM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.1:NumberOfIterations*0.1] );
% grid on

% % for varying VM1
% set(gca,'XTick',[0:4] );
% set(gca,'XTickLabel',[0.1 1 10 100] );
% grid on



subplot(4,1,3)
plot(mct,(final_normoxic_IL),mct,(final_hypoxic_IL),'LineWidth', 2)
% h_legend = legend('1st component','2nd component','Total');
% set(h_legend,'FontSize',14)
xlabel('MCT4_2','FontSize',14)
%xlabel('V_{M_4}','FontSize',14)
%xlabel('K_{M_4}','FontSize',14)
%xlabel('V_{M_1}','FontSize',14)

ylabel('St. st. intra lac','FontSize',14)
title('MCT4 and MCT1 influence on intracellular lactate','FontSize',14)
%title('V_{M_4} influence on intracellular lactate','FontSize',14)
%title('K_{M_4} influence on intracellular lactate','FontSize',14)
%title('V_{M_1} influence on intracellular lactate','FontSize',14)

% for MCT4
set(gca,'XTick',[0:20:400] ); %This are going to be the only values affected.
set(gca,'XTickLabel',[0:0.2:4] ); %This is what it's going to appear in those places.
grid on
% % for varying VM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.01:NumberOfIterations*0.01] );
% grid on

% % for varying KM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.1:NumberOfIterations*0.1] );
% grid on

% % for varying VM1
% set(gca,'XTick',[0:4] );
% set(gca,'XTickLabel',[0.1 1 10 100] );
% grid on


subplot(4,1,4)
plot(mct,(final_normoxic_O),mct,(final_hypoxic_O),'LineWidth', 2)
% h_legend = legend('1st component','2nd component','Total');
% set(h_legend,'FontSize',14)
xlabel('MCT4_2','FontSize',14)
%xlabel('V_{M_4}','FontSize',14)
%xlabel('K_{M_4}','FontSize',14)
%xlabel('V_{M_1}','FontSize',14)

ylabel('St. st. oxygen','FontSize',14)

title('MCT4 and MCT1 influence on oxygen','FontSize',14)
%title('V_{M_4} influence on oxygen','FontSize',14)
%title('K_{M_4} influence on oxygen','FontSize',14)
%title('V_{M_1} influence on oxygen','FontSize',14)

%varying MCT4
set(gca,'XTick',[0:20:400] ); %This are going to be the only values affected.
set(gca,'XTickLabel',[0:0.2:4] ); %This is what it's going to appear in those places.
grid on
% % for varying VM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.01:NumberOfIterations*0.01] );
% grid on

% % for varying VM4
% set(gca,'XTick',[0:10:NumberOfIterations] );
% set(gca,'XTickLabel',[0:10*0.1:NumberOfIterations*0.1] );
% grid on

% % for varying VM1
% set(gca,'XTick',[0:4] );
% set(gca,'XTickLabel',[0.1 1 10 100] );
% grid on











% % to check the stability find jacobian
% syms c1 c2 l1 l2 n1 n2 li1 li2
% 
% jac = jacobian([-omega12*l1 + omega21 * l2 + k1*n1* MCT41*li1 / (KMAX4 + li1) -   k2*n1 * l1 * MCT11  / (KMAX1 + l1) ...
%      - k5 * l1,(-omega21*l2 + omega12 * l1 + k1*n2* MCT42*li2 / (KMAX4 + li2)-  k2*n2 * l2 * MCT12 / (KMAX1 + l2) - k5 * l2),...
%       ((birth(c1) - death(li1,c1)) * n1),((birth(c2) - death(li2,c2)) * n2),S1 - diff12*c1 + diff21*c2 - k6*c1*n1,...
%       S2 + diff12*c1 - diff21 * c2 - k6*c2*n2,(lactate_source1(c1) - k1 * li1 * n1*MCT41/(KMAX4 + li1) + k2*MCT11 * l1*y(1)/(KMAX1 +l1) - k5 * li1),...
%      (lactate_source2(c2) - k1 * li2 * n2*MCT42/(KMAX4 + li2) + k2*MCT11 * l2*y(1)/(KMAX1 +l2) - k5 * li2) ] , [ l1 l2 n1 n2 c1 c2 li1 li2]);
%  % remeber that spaces are quite important here
% 
%  
%  
%  
% jac_ev = subs(jac, [l1 l2 n1 n2 c1 c2 li1 li2], y); % evaluate at the steady state
% jac_ev = double(jac_ev); % make it double 
% % need to check if all eigenvalues are negative (then the poin is stable)
% eigenvalues = (eig(jac_ev));
% 
% j=0;% count number of real eigenvalues
% k=0;% count number of complex eigenvalues
% for i=1:length(eigenvalues)
%     if isreal(eigenvalues(i))==1
%         j=j+1;
%         real_eigen(j) = eigenvalues(i);
%     else
%         k=k+1;
%         complex_eigen(k) = eigenvalues(i);
%     end
% end
% 
% 
% % find minimal eigenvalues
% if j ~= 0
%     max_real = max (real_eigen)
% end
% 
% % if there are complex eigenvalues, check if the real part is negative
% if k ~= 0
%     max_real_part_complex = max(real(complex_eigen))
% end

function value = system(y,MCT41,MCT42, MCT11,MCT12,VM4, KMAX4, VM1,KMAX1)




%    MCT11 = 0.4;
%    MCT12 = 0.4;
%     MCT41 = 0.04;
%     MCT42 = 0.08;
    
    InitialPop =1000;

    k1=1; % this is for intracellulr lactate
    k2=1; % this is for intracellulr lactate
    k5=0.005*InitialPop; % lactate degradation 
    k51=0.005*InitialPop;
    k52=0.5*0.005*InitialPop;

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;

%     omega12=0;
%     omega21=0;

    k4=0.001; % degradation of MCT4
    
    H0=1; % parameter that influences dynamics of MCT4, how fast it increases 
    %as a function of h_i

    %KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    
    % assume it is fixed
    %KMAX4=1; % parameter for lactate dynamics, dependence on MCT4

    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component
    

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    

        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
        h1=exp(beta1*(1-y(5)));
        h2=exp(beta1*(1-y(6)));


    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   

    
    % lactate 1

    value(1) = (-omega12*y(1) + omega21 * y(2) + VM4*y(3)* MCT41 * y(7) / (KMAX4 + y(7)) - ...
        VM1*y(3) * y(1) * MCT11  / (KMAX1 + y(1)) - k5 * y(1)); % assume m1 = 1
    

    % lactate 2      

    value(2) = (-omega21*y(2) + omega12 * y(1) + VM4*y(4)* MCT42 * y(8)/ (KMAX4 + y(8)) - ...
        VM1*y(4) * y(2) * MCT12 / (KMAX1 + y(2)) - k5 * y(2));
        
     
    % population 1

    value(3) = ((birth(y(5)) - death(y(7),y(5))) * y(3));

    % population 2

    value(4) = ((birth(y(6)) - death(y(8),y(6))) * y(4));
    
    % oxygen 1
      
    value(5) = S1 -diff12*y(5) + diff21 * y(6) -k6*y(3)*y(5);
    
    % oxygen 2
    
    value(6) = S2 +diff12*y(5) - diff21 * y(6) -k6*y(4)*y(6);  
    
%     % MCT4 1
% 
%     value(5) = (k3 * h1 /(H0 + h1) - k4 * y(5));
% 
%     % MCT4 2
% 
%     value(6) = (k3 * h2 /(H0 + h2) - k4 * y(6));
%     
%     
%     % intracellular lactate 

   

    value(7) = (lactate_source1(y(5)) - VM4 * MCT41 * y(3)*y(7)/(KMAX4 + y(7)) + VM1*MCT11 * y(3)*y(1)/(KMAX1 +y(1)) - k5 * y(7)); % 1st component
    
    value(8) = (lactate_source2(y(6)) - VM4 * MCT42* y(4)*y(8)/(KMAX4 + y(8)) + VM1*MCT12 * y(4)*y(2)/(KMAX1 +y(2)) - k5 * y(8)); % 1st component

    
   
    
    
    
end

% birth of new cells, dependent on oxygen
function value = birth (oxygen)

    % parameter values from de la Cruz et al. JTB (2016)

    exponent = -0.2;
    a0 = 8.25e3;
    oxycr = 0.02;

    value = a0 * ((oxygen/oxycr)-1)^exponent;
    
    value = 1/value; %%
end
%%%

% death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.2;

    value=nu0*lactate/(L0*oxygen+lactate);
end
%%%



function value = hypoxia (oxygen)

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    

        
        value=exp(beta1*(1-oxygen));
       
end


function value = lactate_source1(oxygen)

            % add source
    S01=10;


    c0=0.05;
        
    % Oxygen dependent source of intracellular lactate (ODIL)
    
    value=S01/(1+(oxygen/c0)); 


end

function value = lactate_source2(oxygen)

            % add source
   
    S02=10;

    c0=0.05;
        
    % Oxygen dependent source of intracellular lactate (ODIL)
    

    value=S02/(1+(oxygen/c0));
end





