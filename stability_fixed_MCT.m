clear all
% stability analysis with varying oxygen but fixed MCT4 and MCT1 different
% Both n1 and n2 are non-zero
% % note that in this case variables are defined two times: outside functions
% and inside


    k4=0.001; % degradation of MCT4
    
    H0=1; % parameter that influences dynamics of MCT4, how fast it increases 
    %as a function of h_i
    
    % parameters for oxygen
    k6 = 1; % how much oxygen is used by cells

    InitialPop = 1000;
    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;
    
    %source of oxygen
    S1=1*InitialPop;
    S2=0;
    % if t<50000
    %   S2=1*InitialPop;
    % else
    %   S2=0;
    % %  S2=1*InitialPop;
    % end
    
    
    

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    
%     if t<75000
%         k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
%         h1=exp(beta1*(1-oxygen1));
%         h2=exp(beta1*(1-oxygen2));

     % parameters which I will use to find steady state lactate
    
    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 
   
   % to check if n=2 is a stable solution
   %stab = birth(oxygen2) - death(l2,oxygen2);
   
   
  
   
   % solve the equation to find steady state oxygen
     
   fun = @oxygen;
   x0 = [2.5,1.5];
    
   [c,fval,exitflag,output]= fsolve(fun,x0)
      

   n1 = pop1(c(1),c(2))
   n2 = pop2(c(1),c(2))
   
   l1 = lac(c(1))
   l2 = lac(c(2))
    
   
   
  %stab = birth(c(1)) - death(l2,c(1))

    % birth of new cells dependent on oxygen
    function value = birth (oxygen)

        % parameter values from de la Cruz et al. JTB (2016)

        exponent = -0.2;
        a0 = 8.25e3;
        oxycr = 0.02;

        value = a0 * ((oxygen/oxycr)-1)^exponent;
        value = 1/value;
    end
    
    
    
%    death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.2;

    value=nu0*lactate/(L0*oxygen+lactate);
end
%


function value = oxygen(c)

    InitialPop = 1000;
    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;


    InitialPop = 1000;
    S1=1*InitialPop;
    k6 = 1;
    
    value(1) = (S1 + diff21 *c(2)-diff12*c(1)) -(k6*c(1))* pop1(c(1),c(2));    
    value(2) = (-diff21*c(2)+diff12*c(1)) -(k6*c(2))* pop2(c(1),c(2)); %first for population 2
    
end

function val = pop1(c1,c2)

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
    InitialPop = 1000;
    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 
    
    mct11 = 0.1;
    mct41 = 2;

    
     val = (-omega12 * lac(c1) + omega21 *lac(c2) - k5 * lac(c1) )/ (-mct41/(KMAX4 + lac(c1)) + lac(c1)*mct11/(KMAX1 +lac(c1)));  
end

function val = pop2(c1,c2)

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
    InitialPop = 1000;
    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 
    
    mct12 = 0.3;
    mct42 = 2;
     

     val = (omega12 * lac(c1) - omega21 *lac(c2) - k5 * lac(c2) )/ (-mct42/(KMAX4 + lac(c2)) + lac(c2)*mct12/(KMAX1 +lac(c2)));  
end



function val = lac(c)

     nu0 = 5e-4;
     L0 = 0.2;

     val = birth(c)*L0*c/(nu0 - birth(c)); 
end


