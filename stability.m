% stability analysis
% First I assumed that there is a non-zero solution
% However, I get wrong solutions even for lactate, so I will check what
% happens when populations die out what is lactate stable solution
% The second part is for n_2 =0 n_1 non-zero
clear all
    InitialPop = 1000;

    oxygen1 = 0.5;
    oxygen2 = 0.5;

    k4=0.001; % degradation of MCT4
    
    H0=1; % parameter that influences dynamics of MCT4, how fast it increases 
    %as a function of h_i

    beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha    
    k3=0.001; % MCT4 dynamics, factor that shows the dependence of hypoxia
    h1=exp(beta1*(1-oxygen1));
    h2=exp(beta1*(1-oxygen2));

    
    % steady state MCT4 value

    mct41 = h1/(H0+h1) * k3/k4
    mct42 = h2/(H0+h2) *k3/k4
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %first part
%     
%     % parameters for lactate dynamics
%     
%     nu0 = 5e-4;
%     L0 = 0.2;
%     
%     % steady state lactate
%     
%     l1 = L0 * oxygen1 * birth(oxygen1)/(nu0 - birth(oxygen1));
%     l2 = L0 * oxygen2 * birth(oxygen2)/(nu0 - birth(oxygen2));
%     birt = birth(oxygen1);
%     
%     % parameters which I will use to find steady state lactate
%     
%     BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
% 
%     % lactate diffusion
%     omega12=InitialPop; 
%     omega21=BiasLactateTransport*InitialPop;
%     
%     KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
%     KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
%     
%     k5=0.005*InitialPop; % lactate degradation 
%     
%     % note that V_M and m1 are all 1
%     
%     
%     n1 = (-omega12 * l1 + omega21*l2 - k5* l1)/ (-mct41/(KMAX4 + l1) + l1/(KMAX1 + l1));
%     
%     n2 = (-omega21 * l2 + omega12*l1 - k5* l2)/ (-mct42/(KMAX4 + l2) + l2/(KMAX1 + l2));
%     
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %second part n_2 zero, n_1 non-zero
 
    nu0 = 5e-4;
    L0 = 0.2;
    
    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation

    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    KMAX1=0.1; % parameter for lactate dynamics, dependence on MCT1
    KMAX4=0.1; % parameter for lactate dynamics, dependence on MCT4
    
    k5=0.005*InitialPop; % lactate degradation 
    
    
    % steady state lactate
    i=0;
    range =0.1:0.05:0.8;
    for oxygen2 = range
        i=i+1;
        l1 = L0 * oxygen1 * birth(oxygen1)/(nu0 - birth(oxygen1))
        l2 = omega12*l1/(omega21 + k5)
        
        % steady state populations
    
        n2 = 0
        n1 = (-omega12*l1 + omega21*l2-k5*l1)/(-mct41/(KMAX4+l1) + l1/(KMAX1 + l1)) 

    
        % to check if n2=0 is a stable solution
        stab(i) = birth(oxygen2) - death(l2,oxygen2);
    end
    
    figure;
    plot(range,stab,'LineWidth',2)
    hold on
    plot([0.5, 0.5],[-2.5*10^-4,0],'r','LineWidth',2)
    plot([0.1 0.5],[0 0],'r','LineWidth',2)
    xlabel('c_2','FontSize',14)
    ylabel('b(c_2) - d(l_2,c_2)','FontSize',14)
    title('Stability of the steady state n_2 = 0','FontSize',14)
    grid on
    
        % birth of new cells dependent on oxygen
    function value = birth (oxygen)

        % parameter values from de la Cruz et al. JTB (2016)

        exponent = -0.2;
        a0 = 8.25e3;
        oxycr = 0.02;

        value = a0 * ((oxygen/oxycr)-1)^exponent;
        value = 1/value;
    end
    
    
    
    % death of cells depending on lactate and oxygen 
function value = death(lactate,oxygen)

    nu0 = 5e-4;
    L0 = 0.2;

    value=nu0*lactate/(L0*oxygen+lactate);
end
%%%
    