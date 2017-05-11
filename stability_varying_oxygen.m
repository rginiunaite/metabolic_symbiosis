clear all
% stability analysis with varying oxygen
% First I will assume n1 non-zero, n2 is zero and try to see if it is
% stable ( n2 = 0 is unstable stab > 0 )
% note that in this case variables are defined two times: outside functions
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
        k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
        h1=exp(beta1*(1-oxygen1));
        h2=exp(beta1*(1-oxygen2));

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
   x0 = 0.5;
   
   c1 = fsolve(fun,x0)
  
   c2 = (diff12/diff21) * c1
   
   n1 = S1/(k6*c1)
   
   
   l1 = lac1(c1)
   l2 = lac2(c1)
    
   
   stab =  (birth(c2) - death(l2,c2))

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


function value = oxygen(x)
    InitialPop = 1000;
    S1=1*InitialPop;
    k6 = 1;
    value = k6 * x * pop(x) -S1;   
end

function val = pop(c)

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
    InitialPop = 1000;
    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 
   

     val = (-omega12 * lac1(c) + omega21 *lac2(c) - k5 * lac1(c) )/ (-1/(KMAX4 + lac1(c))*mct42(c) + lac1(c)/(KMAX1 +lac1(c)));  
end

function val = lac1(c)

     nu0 = 5e-4;
     L0 = 0.2;

     val = birth(c)*L0*c/(nu0 - birth(c)); 
end


function val = lac2(c)

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    InitialPop = 1000;
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 

    nu0 = 5e-4;
    L0 = 0.2;
    
     val = (omega12*birth(c)*L0*c)/((omega21 + k5)*(nu0 - birth(c))); 
end

function val = mct42(c)
k4=0.001;
k3=0.001;
beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha   
H0=1;
     h2 = exp(beta1*(1-c));
     val = h2/(H0+h2) *k3/k4;
end
