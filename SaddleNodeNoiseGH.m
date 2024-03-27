
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code produces one path of the population process x and the solution
%   of the deterministic equation for the strong Allee effect.
%
%   It also generates the corresponding path of the Ornstein-Uhlenbeck
%   process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 1; K = 1;
trho = 0.2501;      % the centered value for the OU process rho (rho star)
lambda = 0.01; phi = 0.00005;     % OU parameters

T = 800; N = 2^15; dt = T/N;

X0 = 0.5; 
rho0 = 0.2501;      % initial value of the process rho

dW = sqrt(dt)*randn(1,N);
W = cumsum(dW);

R = 4; Dt = R*dt; L = N/R;
X = zeros(1,L);
Xtemp = X0;

for j = 1:L
    Winc = sum(dW(R*(j-1)+1:R*j));
    Xtemp = Xtemp + (a*Xtemp^2*(1-Xtemp/K)-(trho+(rho0-trho)*exp(-lambda*j*Dt))*Xtemp)*Dt+(phi/sqrt(2*lambda))*sqrt(1-exp(-2*lambda*j*Dt))*Xtemp*Winc;
    X(j) = Xtemp;
end

rho = zeros(1,L);
theta = zeros(1,L);

integrand = exp(lambda*Dt*(0:L-1));

Winc = zeros(1,L);
for j=1:L    
    Winc(j) = sum(dW(R*(j-1)+1:R*j));
end

for j = 1:L
    rho(j) = trho + (rho0-trho)*exp(-lambda*j*Dt) + phi*exp(-lambda*Dt*j)*sum(integrand(1:j).*Winc(1:j));
    theta(j) = trho+(rho0-trho)*exp(-lambda*Dt*j);
end

rho = [rho0,rho];
theta = [rho0,theta];

Y = zeros(1,L);
Ytemp = X0;
Z = zeros(1,L);
Ztemp = X0;
for j = 1:L
    Ytemp = Ytemp + Dt*(a*Ytemp^2*(1-Ytemp/K)-rho(j)*Ytemp);
    Y(j) = Ytemp;
    Ztemp = Ztemp + Dt*(a*Ztemp^2*(1-Ztemp/K)-rho0*Ztemp);
    Z(j) = Ztemp;
end


EqStable2 = K*0.5*(1+sqrt(1-4*rho/(a*K))); EqUnstable2 = K*0.5*(1-sqrt(1-4*rho/(a*K)));
EqStable2 = EqStable2.*(imag(EqStable2) == 0);
EqUnstable2 = EqUnstable2.*(imag(EqUnstable2) == 0);

clf

subplot(2,1,1) 
hold on
plot(0:Dt:T,EqStable2,'Color',[0.9290 0.6940 0.1250]')
plot(0:Dt:T,EqUnstable2,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')

plot(0:Dt:T,[X0,Y],'r','LineWidth',1.5)
plot(0:Dt:T,[X0,Z],'k','LineWidth',1.5)

xlabel('time')
ylabel('Population density, x')
title('One path of the process x (red) and the solution of the deterministic equation (black)')

subplot(2,1,2)
plot(0:Dt:T,rho)
yline(0.25,'r-')
yline(trho,'-','\rho_{center}')
xlabel('time')
ylabel('\rho')
title('The Ornstein-Uhlenbeck process, \rho')

