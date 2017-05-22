% Variable oxygen model

% Intra- and extra-cellular lactate

% Implicit ODE solver

function msodemodel


IC(1)=1000;
IC(2)=0;
IC(3)=1000;
IC(4)=0;
IC(5)=0.1;
IC(6)=0.1;
%IC(7)=0.7;
%IC(8)=0.3;
% IC(9)=0;
% IC(10)=0;


InitialPop=IC(1);

t0 = 0;
yp0 = [0 0 0 0 0 0  ];
options=odeset('RelTol',1e-6);
% choose suitable initial conditions
[y0,yp0] = decic(@dynsys4,t0,IC,[1 1 1 1 1 1  ],yp0,[0 0 0 0 0 0  ],options,InitialPop);

T=100000;
tspan=[0 T];
y0=IC;
options=odeset('RelTol',1e-4);

[t,y]=ode15i(@dynsys4,tspan,y0,yp0,options,InitialPop);

%figure(4)
%title('Populations')
%plot(t(:),y(:,1),'-g');
%hold on;
%plot(t(:),y(:,3),'-r');
%plot(t(:),y(:,1)+y(:,3),'-k');

%figure(5)
%title('Lactate')
%plot(t(:),y(:,2),'-b');
%hold on;
%plot(t(:),y(:,4),'-m');

%figure(6)
%title('MCT4')
%plot(t(:),y(:,5),'-b');
%hold on;
%plot(t(:),y(:,6),'-m');

%figure(7)
%title('Oxygen')
%plot(t(:),y(:,7),'-b');
%hold on;
%plot(t(:),y(:,8),'-m');

%ages=zeros(1,size(t,1));
%birthrate1=zeros(1,size(t,1));
%birthrate2=zeros(1,size(t,1));

%for i=1:size(t,1)

%  ages(i)=ageg1s(y(i,7));
%  birthrate1(i)=1/(ages(i));
%  ages(i)=ageg1s(y(i,8));
%  birthrate2(i)=1/(ages(i));

%end 

%figure(8)
%plot(t,birthrate1,'-b');
%hold on;
%plot(t,birthrate2,'-m');

figure(11)
subplot(2,2,[1 2])
plot(t(:),y(:,1),'-b');
hold on;
plot(t(:),y(:,3),'-r');
plot(t(:),y(:,1)+y(:,3),'-k');
title('Population');
subplot(2,2,3)
plot(t(:),y(:,2),'-b');
hold on;
plot(t(:),y(:,4),'-r');
title('Lactate');
subplot(2,2,4)
plot(t(:),y(:,5),'-b');
hold on;
plot(t(:),y(:,6),'-r');
title('MCT4');

figure(18)
subplot(3,2,[1 2])
plot(t(:),y(:,1),'-b');
hold on;
plot(t(:),y(:,3),'-r');
plot(t(:),y(:,1)+y(:,3),'-k');
title('Population');
subplot(3,2,3)
plot(t(:),y(:,2),'-b');
hold on;
plot(t(:),y(:,4),'-r');
title('Extracellular lactate');
% subplot(3,2,4)
% plot(t(:),y(:,9),'-b');
% hold on;
% plot(t(:),y(:,10),'-r');
title('Intracellular lactate');
subplot(3,2,5)
plot(t(:),y(:,5),'-b');
hold on;
plot(t(:),y(:,6),'-r');
title('MCT4');
subplot(3,2,6)
% plot(t(:),y(:,7),'-b');
% hold on;
% plot(t(:),y(:,8),'-r');
% title('Oxygen');

% figure(17)
% plot(t(:),y(:,9),'-b');
% hold on;
% plot(t(:),y(:,10),'-r');
% title('Intracellular lactate');

function value=dynsys4(t,y,yp,InitialPop)

oxygen1=0.5;
oxygen2=0.5;


value=zeros(6,1);

k1=1;
k2=1;
k5=0.005*InitialPop;
k51=0.005*InitialPop;
k52=0.5*0.005*InitialPop;

BiasLactateTransport=1;

omega12=InitialPop;
omega21=BiasLactateTransport*InitialPop;

% omega12=0;
% omega21=0;

k4=0.001;
H0=1;

KMAX1=0.1;
KMAX4=0.1;

S1=1*InitialPop;
S2=0;
% if t<50000
%   S2=1*InitialPop;
% else
%   S2=0;
% %  S2=1*InitialPop;
% end

beta1=2.5;
if t<75000
  k3=0.001;
  h1=exp(beta1*(1-oxygen1));
  h2=exp(beta1*(1-oxygen2));
else  
  k3=0.001;  
  %k3=0;
  h1=exp(beta1*(1-oxygen1));
  h2=exp(beta1*(1-oxygen2));
end


k6=1;

BiasOxygenTransport=1;

diffusion12=0.5*InitialPop;
diffusion21=BiasOxygenTransport*0.5*InitialPop;

Therapy=0;
nuTher=1e-4;

%population 1
value(1)=yp(1)-((1/ageg1s(oxygen1)-dr(y(2),oxygen2))*y(1));

%extracellular lactate 1
value(2)=yp(2)-(k1*y(5)*y(1)/(KMAX4+y(2))-k2*y(2)*y(1)/(KMAX1+y(2))-omega12*y(2)+omega21*y(4)-k5*y(2));

%MCT4 1
value(5)=yp(5)-(k3*h1/(H0+h1)-k4*y(5));
%value(5)=0.5*k3-k4*y(5);;

%population 2
value(3)=yp(3)-((1/ageg1s(oxygen2)-dr(y(4),oxygen2))*y(3));

%extracellular lactate 2
value(4)=yp(4)-(k1*y(6)*y(3)/(KMAX4+y(4))-k2*y(4)*y(3)/(KMAX1+y(4))-omega21*y(4)+omega12*y(2)-k51*y(4));

%MCT4 1
value(6)=yp(6)-(k3*h2/(H0+h2)-k4*y(6));
%value(6)=0.5*k3-k4*y(6);

% %Oxygen 1
% value(7)=yp(7)-(S1-k6*y(1)*y(7)-diffusion12*y(7)+diffusion21*y(8));
% 
% %Oxygen 2
% value(8)=yp(8)-(S2-k6*y(3)*y(8)+diffusion12*y(7)-diffusion21*y(8));

%Intracellular lactate 1

% S01=10;
% S02=10;
% 
% c0=0.05;
% 
% Sl1=S01/(1+(y(7)/c0));
% 
% Sl2=S02/(1+(y(8)/c0));
% 
% value(9)=yp(9)-(Sl1-k1*y(5)*y(1)*y(9)/(KMAX4+y(9))+k2*y(2)*y(1)/(KMAX1+y(2))-k5*y(9));
% 
% value(10)=yp(10)-(Sl2-k1*y(6)*y(3)*y(10)/(KMAX4+y(10))+k2*y(4)*y(3)/(KMAX1+y(4))-k52*y(10));

%%%

function value=ageg1s(oxygen)

% parameter values from de la Cruz et al. JTB (2016)

exponent = -0.2;
a0 = 8.25e3;
oxycr=0.02;

value = a0*((oxygen/oxycr)-1)^exponent;

%%%

function value=dr(lactate,oxygen)

nu0=5e-4;
L0=0.1;

value=nu0*lactate/(L0*oxygen+lactate);

%%%

function value=drther(lactate,oxygen)

nu0=5e-4;
L0=0.2;

Therapy=0;
nuther=1e-4;

value=Therapy*nuther+nu0*lactate/(L0*oxygen+lactate);
