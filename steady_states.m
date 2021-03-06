% this is to find steady states numerically, all before was wrong
% a little bit approximate value
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
    % fixed oxygen (temporal)
    
    oxygen1= 0.5;
    oxygen2 =0.5;


 MCT11 = 1;
    MCT12 = 1;
    
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

    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4

    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component
    

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    

        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia



    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









IC(1) = 0.1;        % Initial extracellular lactate, 1st compartment
IC(2) = 0.1;        % Initial extracellular lactate, 2nd compartment
IC(3) = 1000;     % Initial population 1st compartment
IC(4) = 1000;     % Initial population 2nd compartment
IC(5) = 0.7;      % Initial oxygen 1st compartment    
IC(6) = 0.3;      % Initial oxygen 2nd compartment  
IC(7) = 0.1;      % Initial MCT4 1st comparmtent
IC(8) = 0.1;      % Initial MCT4 2nd compartment
%IC(9)=0;
%IC(10)=0;
InitialPop = IC(3); % have variable for initial population in the first compartment

i=0;
for k3 = [0.001,0.005,0.01,0.05]
i=i+1;
   fun = @(y) system(y,k3);
  
     options=optimset('disp','iter','LargeScale','off','TolFun',1e-10,'MaxIter',1e8,'MaxFunEvals',1e8,'TolX',1e-10);
 
    [y,fval,exitflag,output] = fsolve(fun,IC,options);
    
len = length(y(:,3)); % number of simulations done

final_normoxic(i) = y(len,3);
final_hypoxic(i) = y(len,4);
final_total(i) = y(len,3) + y(len,4);
end
%system(IC)

range=1:i;

figure

% Population
subplot(2,1,1);
plot(range,(final_normoxic),'LineWidth', 2)
xlabel('k_3','FontSize',14)
ylabel('Normoxic population','FontSize',14)
set(gca,'XTick',[1 2 3 4] ); %This are going to be the only values affected.
set(gca,'XTickLabel',[0.001 0.005 0.01 0.05] ); %This is what it's going to appear in those places.
grid on


subplot(2,1,2)

plot(range,final_hypoxic,'LineWidth', 2)
xlabel('k_3','FontSize',14)
ylabel('Hypoxic population','FontSize',14)
%h_legend = legend('1st component','2nd component','Total')
%set(h_legend,'FontSize',14)

%title(['k_3 influence on total population'],'FontSize',14)
set(gca,'XTick',[1 2 3 4] ); %This are going to be the only values affected.
set(gca,'XTickLabel',[0.001 0.005 0.01 0.05] ); %This is what it's going to appear in those places.
grid on








% to check the stability find jacobian
syms l1 l2 n1 n2 c1 c2 m41 m42

jac = jacobian([lactate_source1(c1)-omega12*l1 + omega21 * l2 + k1*n1* m41 / (KMAX4 + l1) -   k2*n1 * l1 * MCT11  / (KMAX1 + l1) ...
     - k5 * l1,(lactate_source2(c2)-omega21*l2 + omega12 * l1 + k1*n2* m42 / (KMAX4 + l2)-  k2*n2 * l2 * MCT12 / (KMAX1 + l2) - k5 * l2),...
      ((birth(c1) - death(l1,c1)) * n1),((birth(c2) - death(l2,c2)) * n2),S1-diff12*c1 + diff21 * c2-k6*n1*c1,S2+diff12*c1 - diff21 * c2-k6*n2*c2,...
      (k3 * hypoxia(c1) /(H0 + hypoxia(c1)) - k4 * m41),...
      (k3 * hypoxia(c2) /(H0 + hypoxia(c2)) - k4 * m42) ] , [l1 l2 n1 n2 c1 c2 m41 m42]);


M = subs(jac, [l1 l2 n1 n2 c1 c2 m41 m42], y); % evaluate at the steady state

% need to check if all eigenvalues are negative (then the poin is stable)
eigenvalues = (eig(M));
M = double(M); % make it double 
j=0;% count number of real eigenvalues
k=0;% count number of complex eigenvalues
for i=1:length(eigenvalues)
    if isreal(eigenvalues(i))==1
        j=j+1;
        real_eigen(j) = eigenvalues(i);
    else
        k=k+1;
        complex_eigen(k) = eigenvalues(i);
    end
end


% find minimal eigenvalues
max_eigen = max (real_eigen)
% if there are complex eigenvalues, check if the real part is negative
if k ~= 0
    max_real_part_complex = max(real(complex_eigen))
end

function value = system(y,k3)

    MCT11 = 1;
    MCT12 = 1;
    
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

    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4

    S1=1*InitialPop; % source of oxygen depends on how much population I have initially 
    S2 = 0; % when anitangiogenesis treatment is applied, no oxygen source to the second component
    

    

    %    k3=0.01; % MCT4 dynamics, factor that shows the dependence of hypoxia



    k6 = 1; % how much oxygen is used by cells

    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
   

    
    % lactate 1

    value(1) = (lactate_source1(y(5))-omega12*y(1) + omega21 * y(2) + k1*y(3)* y(7) / (KMAX4 + y(1)) - ...
        k2*y(3) * y(1) * MCT11  / (KMAX1 + y(1)) - k5 * y(1)); % assume m1 = 1
    

    % lactate 2      

    value(2) = (lactate_source2(y(6))-omega21*y(2) + omega12 * y(1) + k1*y(4)* y(8) / (KMAX4 + y(2)) - ...
        k2*y(4) * y(2) * MCT12 / (KMAX1 + y(2)) - k5 * y(2));
        
     
    % population 1

    value(3) = ((birth(y(5)) - death(y(1),y(5))) * y(3));

    % population 2

    value(4) = ((birth(y(6)) - death(y(2),y(6))) * y(4));
    
    % oxygen 1
      
    value(5) = S1 -diff12*y(5) + diff21 * y(6) -k6*y(3)*y(5);
    
    % oxygen 2
    
    value(6) = S2 +diff12*y(5) - diff21 * y(6) -k6*y(4)*y(6);  
    
    % MCT4 1

    value(7) = (k3 * hypoxia(y(5)) /(H0 + hypoxia(y(5))) - k4 * y(7));

    % MCT4 2

    value(8) = (k3 * hypoxia(y(6)) /(H0 + hypoxia(y(6))) - k4 * y(8));
    
   
    
    
    
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
    S01=100;

    c0=0.05;
        
    % Oxygen dependent source of intracellular lactate (ODIL)
    
    value=S01/(1+(oxygen/c0)); 


end

function value = lactate_source2(oxygen)

    % add source
   
    S02=100;

    c0=0.05;
        
    % Oxygen dependent source of intracellular lactate (ODIL)
    

    value=S02/(1+(oxygen/c0));
end
