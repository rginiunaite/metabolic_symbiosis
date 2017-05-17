clear all
% stability analysis with intracellular lactate
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
    S2=0;%0.01*InitialPop;
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
   i = 0; % for counting
   range = 0:5:150;
   %for S2 = range  
   S2 = 0;
       i=i+1;
   fun = @(c) oxygen(c,S2);
   x0 = [0.7,0.3];
    
   [c,fval,exitflag,output]= fsolve(fun,x0)
      
   n1(i) = pop1(c(1),c(2))
   n2(i) = pop2(c(1),c(2))

   
   l_intra1 = lac_intra(c(1))
   l_intra2 = lac_intra(c(2))
    
   mct41 = mct4((c(1)))
   mct42 = mct4((c(2)))
   
   

   %end
   
%    plot(range,log(n1),range,log(n2),'Linewidth',2)
%    h_legend = legend('Population 1','Population 2')
%    set(h_legend,'Fontsize',14)
%    xlabel('S2','Fontsize',14)
%    ylabel('log of steady state populations','Fontsize',14)
   
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


function value = oxygen(c,S2)

    InitialPop = 1000;
    BiasOxygenTransport=1;

    diff12 = 0.5*InitialPop; % diffusion depends on Initial Poulation in the 1st comp
    
    diff21 = BiasOxygenTransport*0.5*InitialPop;


    InitialPop = 1000;
    S1=1*InitialPop;
    k6 = 1;
    %S2 = 0.001*InitialPop;
    %S2 = 150;
    
    value(1) = (S1 + diff21 *c(2)-diff12*c(1)) -(k6*c(1))* pop1(c(1),c(2));    
    value(2) = (S2 -diff21*c(2)+diff12*c(1)) -(k6*c(2))* pop2(c(1),c(2)); %first for population 2
    
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
    k51=0.005*InitialPop;
    
     val = (source1(c1)- k51 * lac_intra(c1))/(mct4(c1)*lac_intra(c1)/(KMAX4 + lac_intra(c1)) - lac_extra1(c1,c2)/(KMAX1 + lac_extra1(c1,c2)) );
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
    k51=0.005*InitialPop;

     val = (source2(c2)- k51 * lac_intra(c2))/(mct4(c2)*lac_intra(c2)/(KMAX4 + lac_intra(c2)) - lac_extra2(c1,c2)/(KMAX1 + lac_extra2(c1,c2)) );
end

function val = lac_intra(c)

     nu0 = 5e-4;
     L0 = 0.2;

     val = birth(c)*L0*c/(nu0 - birth(c)); 
end

function val = lac_extra1(c1,c2)

    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
    InitialPop = 1000;
    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    k5=0.005*InitialPop; % lactate degradation 
    
    k51=0.005*InitialPop;


    val = (omega21*(-source1(c1) + k51 *lac_intra(c1) - source2(c2) + k51 * ...
        lac_intra(c2)) +k5*(-source1(c1) +k51*lac_intra(c1)))/(-k5*omega21 - k5 * omega12 - k5^2);

end

function val = lac_extra2(c1,c2)


    BiasLactateTransport=1; % initially assume there is no bias in lactate transportation
    InitialPop = 1000;
    % lactate diffusion
    omega12=InitialPop; 
    omega21=BiasLactateTransport*InitialPop;
    
    k5=0.005*InitialPop; % lactate degradation 
    
    k51=0.005*InitialPop;

    val = (omega12 * lac_extra1(c1,c2) + k5 * lac_extra1(c1,c2) - source1(c1) + k51 * lac_intra(c1))/omega21;

end


function val = mct4(c)

k4=0.001;
k3=0.001;
beta1 = 2.5; %parameter for the hypoxia dependent on HIF-1alpha   
H0=1;
     h2 = exp(beta1*(1-c));
     val = h2/(H0+h2) *k3/k4;
     
end

    % Oxygen dependent source of intracellular lactate (ODIL)
    
    function  val = source1(c)
    
        S01=10;
        c0=0.05;
        val = S01/(1+c/c0);
    
    end

    function  val = source2(c)
    
        S02=10;
        c0=0.05;
        val = S02/(1+c/c0);
    
    end




