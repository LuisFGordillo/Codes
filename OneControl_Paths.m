clearvars
warning off

a = 1; K = 1;
trho = 0.247;      % the centered value for the OU process rho (rho star)
lambda = 0.01; 
phi = 0.04;     % OU parameters

psi = 0.1; 
eps = 0.005;

P = 0.95;
k =-log(1-P);
T = 40000; N = 2^17; dt = T/N;
X0 = 0.5*(1+sqrt(1-4*trho)); 
rho0 = trho;                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial value of the process rho


xs = 0.5*(1+sqrt(1-4*trho));                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stable equilibrium
xu = 0.5*(1-sqrt(1-4*trho));                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stable equilibrium
J = [2*xs-3*xs^2-trho -xs; 0 -lambda];      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Jacobian aevaluated at xs
S = [psi*xs 0; 0 phi^2];                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix sigma*sigma'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Increments of Brownian motion
dW = sqrt(dt)*randn(2,N);
R = 4; Dt = R*dt; L = N/R;



rho = zeros(1,L);                           
rhotemp = rho0;
X = zeros(1,L);                             
Xtemp = X0;
 
% Get the dimensions of matrix J
 n = size(J, 1);
 
 % This equation solves vec(J * W + W * J') = vec(-S)
 I = eye(n);
 lhs = kron(I,J) + kron(J,I);  % Left-hand side matrix for the linear system
 rhs = -S(:);                   % Right-hand side vectorized
 
 % Solve the linear system
 wtemp = lhs \ rhs;
 
 % Reshape the solution vector back into the matrix form
 SM = reshape(wtemp, n, n);
 
 SM(1,1) = 0.8;
 SM(2,2) = SM(1,1);



B = [1;0];
BB = [1 0];
P = eye(2)-B*BB;


CM = -0.5*BB*(S+J*SM+SM*J')*(P+eye(2))*inv(SM);

Winc = zeros(2,L);
for j=1:L    
    Winc(1,j) = sum(dW(1,R*(j-1)+1:R*j));
    Winc(2,j) = sum(dW(2,R*(j-1)+1:R*j));
end

for j = 1:L    
    Xtemp = Xtemp + Dt*(a*Xtemp^2*(1-Xtemp/K)-rhotemp*Xtemp)+eps*sqrt(psi*Xtemp)*Winc(1,j) + Dt*(CM(1,1)*(Xtemp-xs) + CM(1,2)*(rhotemp-trho));
    rhotemp = rhotemp - lambda*(rhotemp - trho)*Dt + eps*phi*Winc(2,j);% + Dt*(CM(2,1)*(Xtemp-xs) + CM(2,2)*(rhotemp-trho));

    X(j) = Xtemp;
    rho(j) = rhotemp;
end

rho = [rho0,rho];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create plots of paths
figure(3)
hold on
X=[X0,X];
X2D = real(X);

plot(rho,X,'MarkerEdgeColor',"#A2142F")
ylabel('x','FontSize',16);
xlabel('\rho','FontSize',16);
title('(b) Dynamics with control','FontSize',16)


%axis([0.2485 0.251 0.44 0.54])
theta=0.240:0.000001:0.25;
EqStable1 = K*0.5*(1+sqrt(1-4*theta/(a*K))); EqUnstable1 = K*0.5*(1-sqrt(1-4*theta/(a*K)));


z = length(rho);
TimeMarks = 150;
rhotemp=rho(1:TimeMarks:z);
Ztemp=X(1:TimeMarks:z);
plot(rhotemp,Ztemp,'g.')

plot(theta,EqStable1,'r-')
plot(theta,EqUnstable1,'r--')
xline(trho,'LineWidth',1.5);
xline(0.25,'--','LineWidth',1.5)
eq = K*0.5*(1+sqrt(1-4*trho/(a*K)));
plot(trho, eq,'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the confidence ellipse

plot(trho, xu,'ko','MarkerFaceColor', 'none', 'MarkerSize', 9)

E = inv(SM);

ellipse = @(ou,pop) E(1,1)*(pop-xs).^2+E(2,2)*(ou-trho).^2+(E(2,1)+E(1,2))*(pop-xs).*(ou-trho)-2*eps^2*k;
fimplicit(ellipse,'r','LineWidth',3)

axis([0.243 0.255 0.43 0.6]);


%%%% Plot the separatrix %%%%

yRange = 0.246:0.0003:0.251;
xRange = 0.4:0.01:0.55;

% Create a meshgrid for x and y
[X, Y] = meshgrid(xRange, yRange);

% Define the vector field functions f(x, y) and g(x, y)
f = @(x, y) x.^2.*(1-x)-y.*x;    % Example function for x'
g = @(x, y) -lambda*(y-trho);    % Example function for d/dt(rho)

% Define the system of differential equations for ode45
dynamics = @(t, z) [f(z(1), z(2)); g(z(1), z(2))];

% Identify the saddle point and its Jacobian
saddle_point = [xu,trho]; % Known saddle point for this example
J = [2*xu-3*xu^2-trho -xu; 0 -lambda];      % Jacobian matrix at the saddle point

% Find the eigenvalues and eigenvectors of the Jacobian
[VV, D] = eig(J);
stable_direction = VV(:,2); % Eigenvector for the unstable eigenvalue

% Perturbation size
epsilon = 0.0005;

% Compute initial points near the saddle point in the stable direction
z0_pos = saddle_point + epsilon * stable_direction';
z0_neg = saddle_point - epsilon * stable_direction';

% Integrate forward and backward from these points to plot the stable separatrix

[T2, Z2] = ode45(dynamics, [0, -6000], z0_pos);   % Backward in time
[T4, Z4] = ode45(dynamics, [0, -1000], z0_neg);   % Backward in time

plot(Z2(:, 2), Z2(:, 1), 'b', 'LineWidth', 1.5); % Positive direction backward
plot(Z4(:, 2), Z4(:, 1), 'b', 'LineWidth', 1.5); % Negative direction backward

hold off;

ax = gca;
ax.YTick = linspace(0.43, 0.6, 5);
ax.XTick = linspace(0.243, 0.2550, 3);

